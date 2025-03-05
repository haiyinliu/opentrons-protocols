from opentrons import protocol_api
from opentrons import types
from opentrons.protocol_api.labware import OutOfTipsError

metadata = {
    "protocolName": "Reverse transcription and strand-switching (SQK-LSK114)",
    "description": "Based on the direct cDNA sequencing protocol (DCS_9187_v114_revJ_12Dec2024)",
    "author": "Haiyin Liu"
    }

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.18",
}

def run(protocol: protocol_api.ProtocolContext):
    #======== PARAMETERS ========
    # SAMPLE PARAMETERS
    NUM_SAMPLES = 8  # Define the number of samples (up to 48)
    NUM_COLUMNS = (NUM_SAMPLES + 7) // 8  # Calculate number of columns for 96-well plates
    ELUTION_VOL = 21    # µl of NFW to resuspend the beads in at the last elution

    # TESTING PARAMETERS
    DRY_RUN = 1      # use 0.01 to shorten wait times if it is dry run, otherwise 1
    BEAD_RUN = 1        # use 0.01 for testing without beads, otherwise 1 

    skip_firststrand = False     # Toggle when testing certain blocks, same below
    skip_65_5min = False
    skip_strandswitch = False
    skip_42_90min = False
    skip_rnase = False
    skip_bead17 = False
    skip_wash1 = False
    skip_elution = False
    skip_2ndstrand = False
    skip_bead40 = False
    skip_wash2 = False
    skip_finalelution = False

    # VOLUME AND DISTANCE SETTINGS
    deadvol_reservoir = 1500
    deadvol_plate = 10
    clearance_reservoir = 2
    clearance_bead_pellet = 1.5
    clearance_beadresuspension = 3
    offset_x = 1

    #======== DECK SETUP ========
    ## COLUMN 1 
    thermocycler = protocol.load_module('thermocycler module gen2')
    
    # Load the polyA plate on the heater shaker module (just as a holding location)
    heatershaker = protocol.load_module('heaterShakerModuleV1','C1')
    hs_adapter = heatershaker.load_adapter('opentrons_96_pcr_adapter')
    heatershaker.open_labware_latch()
    plate1_polyA = hs_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    protocol.pause("Place polyA plate on the Heater-Shaker Module and resume.")
    heatershaker.close_labware_latch()

    temp_block = protocol.load_module('temperature module gen2', 'D1')  
    temp_adapter = temp_block.load_adapter('opentrons_96_well_aluminum_block')     
    cold_plate = temp_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")

    ## COLUMN 2 
    mag_block = protocol.load_module('magneticBlockV1', 'A2')
    plate3_eluate = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "B2")
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    plate2_sample = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "D2")
    
    ## COLUMN 3 
    tiprack_50ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'A3')
    tiprack_50ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B3')
    tiprack_200ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C3')
    waste_chute = protocol.load_waste_chute()
    
    ## COLUMN 4 - staging area
    tiprack_50ul_3 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'A4')
    tiprack_50ul_4 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B4')
    tiprack_200ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C4')
    # D4 is an empty slot to move tip boxes to
 
    #======== PIPETTES ========
    # load active tips for each pipette
    p50m = protocol.load_instrument('flex_8channel_50', 'right')
    p1000m = protocol.load_instrument('flex_8channel_1000', 'left')
    
    # define "active" and "staging" tip boxes for each pipette
    p50m.tip_racks = [tiprack_50ul_1, tiprack_50ul_2]
    p50_staging = [tiprack_50ul_3, tiprack_50ul_4]
    
    p1000m.tip_racks = [tiprack_200ul_1]
    p1000_staging = [tiprack_200ul_2]

    #======== REAGENT WELLS ========
    # define well locations for where samples are moved from/to 
    polyA_wells = plate1_polyA.columns()[6 : 6 + NUM_COLUMNS]  # oligo dT-bead eluates from polyA protocol, old plate columns 7-12)
    firststrand_wells = plate2_sample.columns()[:NUM_COLUMNS]  # wells for the first strand synthesis
    secondstrand_wells = plate2_sample.columns()[6 : 6 + NUM_COLUMNS] # wells for the second strand synthesis
    eluate_wells = plate3_eluate.columns()[:NUM_COLUMNS]  # wells for the final eluate

    # cold plate (for reagents and the eluate for this protocol)
    MM1 = cold_plate.columns_by_name()['1']
    MM2 = cold_plate.columns_by_name()['2']
    MM3 = cold_plate.columns_by_name()['3']
    beads = [cold_plate.columns_by_name()[name] for name in ['4', '5']]
    rnase = cold_plate.columns_by_name()['6']
    
    eluate = [
        cold_plate.columns_by_name()[name] for name in ['7', '8', '9', '10', '11', '12']]

    # beads (2X 47 uL per sample)
    for index, column in enumerate(beads):
        if not index:
            num = NUM_COLUMNS if NUM_COLUMNS <= 3 else 3    # only 1 column for up to 3 samples
        else:
            num = NUM_COLUMNS - 3 if NUM_COLUMNS > 3 else 0 # 2 columns for more than 3 samples
        # sets volume needed in each column
        # note: I don't know where this is used if at all
        if num:
            column[0].liq_vol = 47*num + deadvol_plate
        else:
            column[0].liq_vol = 0

    # reagent reservoir
    ethanol = [reservoir.wells_by_name()[well] for well in ['A1', 'A2', 'A3']]
    for well in ethanol: well.liq_vol = 4*200*2*p1000m.channels + deadvol_reservoir
    NFW = reservoir['A5']

    # waste 
    waste = [reservoir.wells_by_name()[well] for well in ['A9','A10', 'A11', 'A12']]

    #======== DEFINING LIQUIDS ========   
    ethanol_liquid = protocol.define_liquid(
        name="Ethanol",
        description="80% Ethanol in NFW in reservoir",
        display_color="#008080"  # Teal
    )

    MM1_liquid = protocol.define_liquid(
    name="MM1",
    description="VN primer and dNTP mix, 3.5µl per sample",
    display_color="#4682B4"  # SteelBlue
    )

    MM2_liquid = protocol.define_liquid(
        name="MM2",
        description="RT buffer, RNAseOUT, StrSw Primer and Maxima RTase, 8µl per sample",
        display_color="#DAA520"  # GoldenRod
    )

    MM3_liquid = protocol.define_liquid(
        name="MM3",
        description="LongAmp Taq and PR2 primer",
        display_color="#50C878"  # EmeraldGreen
    )

    beads_liquid = protocol.define_liquid(
        name="Beads",
        description="Magnetic beads, 57µl per sample",
        display_color="#99664B"  # AmethystPurple
    )

    NFW_liquid = protocol.define_liquid(
        name="NFW",
        description="Nuclease Free Water, 20µl per sample",
        display_color="#999999"  # Gray
    )

    polyA_samples_liquid = protocol.define_liquid(
        name="PolyA Samples",
        description="PolyA RNA samples",
        display_color="#DC143C"  # Crimson
    )

    eluate_liquid = protocol.define_liquid(
        name="Eluate",
        description="Final eluate",
        display_color="#FFD700"  # Gold
    )

    # Define volumes
    ethanol_vol = 0.2 * 4 * NUM_COLUMNS * 8 + deadvol_reservoir
    MM1_vol = 4.5 * NUM_COLUMNS + deadvol_plate
    MM2_vol = 9 * NUM_COLUMNS + deadvol_plate
    MM3_vol = 30 * NUM_COLUMNS + deadvol_plate
    beads_vol = 57  * NUM_COLUMNS + deadvol_plate
    NFW_vol = 20 * NUM_COLUMNS + deadvol_plate
    polyA_samples_vol = 15
    eluate_vol = 20

    #======== LOADING LIQUIDS ========
    # sample plate
    for column in polyA_wells:
        for well in column:
            well.load_liquid(liquid=polyA_samples_liquid, volume=polyA_samples_vol)

    # reservoir
    for column in reservoir.columns()[0:2]:
        for well in column: 
            well.load_liquid(liquid=ethanol_liquid, volume=ethanol_vol)

    # cold plate
    # dictionary for mapping liquids on the plate 
    liquid_mapping = {
        # Master Mixes in columns "1", "2", and "3"
        '1': (MM1_liquid, MM1_vol),
        '2': (MM2_liquid, MM2_vol),
        '3': (MM3_liquid, MM3_vol),

        # beads goe in columns "4" and "5"
        '4': (beads_liquid, beads_vol),
        '5': (beads_liquid, beads_vol),

        # NFW goes in column "6"
        '6': (NFW_liquid, NFW_vol),

        # Eluate goes in columns "7" through "12"
        '7': (eluate_liquid, eluate_vol),
        '8': (eluate_liquid, eluate_vol),
        '9': (eluate_liquid, eluate_vol),
        '10': (eluate_liquid, eluate_vol),
        '11': (eluate_liquid, eluate_vol),
        '12': (eluate_liquid, eluate_vol)
    }

    # loop over the mapping to load liquids into column
    for col_name, (liquid, vol) in liquid_mapping.items():
        for well in cold_plate.columns_by_name()[col_name]:
            well.load_liquid(liquid=liquid, volume=vol)

    #======== RUN SETUP ========
    # HELPER FUNCTIONS 
    # FUNCTION 1: PIPETTE PICK UP
    def pick_up_or_refill(pipette):
        """Pick up a tip or pause for replacement if needed."""
        try:
            pipette.pick_up_tip()
        except OutOfTipsError:
            pause_attention(
             """Please Refill the {} Tip Boxes
                and Empty the Tip Waste.""".format(current_pipette))
            pipette.reset_tipracks()
            pipette.pick_up_tip()
  
    # FUNCTION 2: SLOW PIPETTE WITHDRAWAL    
    def slow_tip_withdrawal(pipette, well, z=0, delay_seconds=0):
        pipette.default_speed /= 10
        if delay_seconds > 0:
            protocol.delay(seconds=delay_seconds)
        pipette.move_to(well.top(z))
        pipette.default_speed *= 10
    
    # FUNCTION 3: SIDE TOUCH WITH BLOWOUT
    def side_touch_w_blowout(pipette,well,pos=-1):
            pipette.default_speed /= 40
            pipette.move_to(well.top(pos).move(types.Point(x=well.diameter / 2, y=0, z=0))) 
            pipette.blow_out()
            protocol.delay(seconds=1.5)
            pipette.default_speed *= 2
            pipette.move_to(well.top())
            pipette.default_speed *= 20

    # FUNCTION 4: BEAD MIXING (for when beads are in solution)
    def bead_mixing(pipette, reps, vol, speed=0.5):
        for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
            if not pipette.has_tip: 
                pick_up_or_refill(pipette)
            
            # slow mixes
            pipette.mix(reps, vol, column[0].bottom(3), rate=speed)
            slow_tip_withdrawal(pipette, column[0], -2)
            
            # side touch with blowout after last mix
            side_touch_w_blowout(pipette,column[0],pos=-5)                  
            pipette.drop_tip()

    # FUNCTION 5: REMOVE SUPERNATANT
    def remove_sup(pipette,volume1,volume2,waste_well):       
        pick_up_or_refill(pipette)
        
        pipette.move_to(column[0].top())
        pipette.air_gap(20)
        pipette.aspirate(volume1, column[0].bottom(4), rate=0.1)        
        pipette.aspirate(volume2, column[0].bottom(clearance_bead_pellet), rate=0.02)
        protocol.delay(seconds=1)

        pipette.dispense(volume1+volume2+20, waste[waste_well].top(), rate=2)
        protocol.delay(seconds=1)
        pipette.blow_out()

    # FUNCTION 6: WASH MIXING (for when beads are in a pellet)
    def wash_mixing(pipette, volume,reps):
        
        # rotating locations to rinse down the bead "ring" from all sides
        locations = [
            types.Point(x=offset_x, y=0, z=0),
            types.Point(x=-offset_x, y=0, z=0),
            types.Point(x=0, y=offset_x, z=0),
            types.Point(x=0, y=-offset_x, z=0)]

        for point in locations:
            loc = column[0].bottom(3).move(point)
            pipette.mix(3, volume / 2, loc, rate=2)
        
        # slower mixes with whole volume
        pipette.mix(reps, volume, column[0].bottom(3), rate=0.5)
        slow_tip_withdrawal(pipette, column[0], -2)
        
    #======== STEP BY STEP PROTOCOL ========
    # Set the temperature of the cooler block and thermo cycler
    temp_block.set_temperature(celsius=4)  
    thermocycler.set_block_temperature(temperature=20)
    
    ## 4. Reverse transcription and strand-switching
    # STEP 4.3 - mix RNA input and VN primer mix
    if not skip_firststrand:
        
        # distribute MM1 to required wells on the sample plate
        p50m.distribute(
            volume=3.5,
            source=MM1[0].bottom(1),
            dest=[col[0].bottom(1) for col in firststrand_wells[:NUM_COLUMNS]],
            new_tip='once',
            aspirate_rate=0.5,
            dispense_rate=1,
            blow_out=True,
            blow_out_location='destination well'
        )
        
        # Transfer polyA samples to sample plate and mix
        p50m.transfer(
            volume=7.5,
            source=[col[0].bottom(1) for col in polyA_wells],
            dest=[col[0].bottom(1) for col in firststrand_wells[:NUM_COLUMNS]],
            new_tip='always',
            aspirate_rate=0.5,
            dispense_rate=1,
            mix_after=(3, 7.5),  # Mix 3 times with 7.5µL after dispensing
            blow_out=True,
            blow_out_location='destination well',
            touch_tip=12    # Touch tip 12mm below the top of the well, check if this is appropriate
        )

    # STEP 4.5 - heat to 65ºC for 5 minutes and immediately move to cold plate
    if not skip_65_5min:
        thermocycler.set_lid_temperature(80)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(65, hold_time_minutes=5*DRY_RUN, block_max_volume=15)
        thermocycler.set_block_temperature(4, hold_time_minutes=1*DRY_RUN)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)
        thermocycler.set_block_temperature(20)
            
    # STEP 4.6-4.11 - add strand switch primer mix and Maxima H Minus Rev transcriptase
    if not skip_strandswitch:
        p50m.transfer(
            volume=9,
            source=MM2[0].bottom(1),
            dest=[col[0].bottom(1) for col in firststrand_wells[:NUM_COLUMNS]],
            new_tip='always',
            aspirate_rate=0.5,
            dispense_rate=1,
            mix_after=(3, 10),  # Mix 3 times with 7.5µL after dispensing
            blow_out=True,
            blow_out_location='destination well',
            touch_tip=12    # Touch tip 12mm below the top of the well, check if this is appropriate
        )

    # STEP 4.12 - incubate at 42ºC for 90 minutes and 85ºC for 5 minutes
    if not skip_42_90min:
        thermocycler.set_lid_temperature(95)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(42, hold_time_minutes=90*DRY_RUN, block_max_volume=20)
        thermocycler.set_block_temperature(85, hold_time_minutes=5*DRY_RUN, block_max_volume=20)
        thermocycler.set_block_temperature(4)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    ## 5. RNA degradation and second strand synthesis
    # STEP 5.3-5.4 - add RNase and incubate at 37ºC for 10 minutes
    if not skip_rnase:
        # Transfer rnase to clean wells and mix
        p50m.transfer(
            volume=1,
            source=rnase[0].bottom(1),
            dest=[col[0].bottom(1) for col in firststrand_wells[:NUM_COLUMNS]],
            new_tip='always',
            aspirate_rate=0.5,
            dispense_rate=1,
            mix_after=(2, 10),  # Mix 3 times with 10µL after dispensing
            mix_kwargs={'aspirate_rate': 0.5, 'dispense_rate': 1},
            blow_out=True,
            blow_out_location='destination well',
            touch_tip=10    # Touch tip 12mm below the top of the well, check if this is appropriate
        )

        thermocycler.set_lid_temperature(55)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(37, hold_time_minutes=10*DRY_RUN, block_max_volume=20)
        thermocycler.set_block_temperature(20)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    # STEP 5.5+5.7 - Mix beads before adding to the plate
    if not skip_bead17:
        # Calculate mixing volumes for each bead column based on NUM_COLUMNS
        bead_mix_vols = []
        if NUM_COLUMNS <= 3:
            # Only use first column of beads
            bead_mix_vols.append(50 * NUM_COLUMNS)  # 50µl per sample column
            bead_mix_vols.append(0)  # Second column not used
        else:
            # Use both bead columns
            bead_mix_vols.append(150)  # First column always 150µl for 3 columns
            bead_mix_vols.append(50 * (NUM_COLUMNS - 3))  # Second column for remaining columns

        # Resuspend beads in each column with p1000
        for idx, (bead_col, mix_vol) in enumerate(zip(beads, bead_mix_vols)):
            if mix_vol > 0:  # Only process columns with beads
                p1000m.pick_up_tip()
                
                # Thorough mixing at different heights
                p1000m.mix(5, mix_vol * 0.8, bead_col[0].bottom(2))  # Bottom mixing
                p1000m.mix(5, mix_vol * 0.8, bead_col[0].bottom(5))  # Middle mixing
                
                # Final slow mixes
                for _ in range(3):
                    p1000m.aspirate(mix_vol, bead_col[0].bottom(2), rate=0.5)
                    p1000m.dispense(mix_vol, bead_col[0].bottom(5), rate=0.5)
                
                # Slow withdrawal
                slow_tip_withdrawal(p1000m, bead_col[0], -2)
                p1000m.drop_tip()

        # Add 17µl of beads to sample columns (with pre- and post-mixing)
        for col_idx in range(NUM_COLUMNS):
            bead_source = beads[0] if col_idx < 3 else beads[1]  # Switch source after 3 columns
            
            p50m.pick_up_tip()
            p50m.mix(3, 50, bead_source[0].bottom(2), rate=0.5)  # Gentle remix before aspiration
            
            # Transfer beads with mixing
            p50m.transfer(
                17,
                bead_source[0].bottom(2),
                firststrand_wells[col_idx][0].bottom(2),
                new_tip='never',
                mix_after=(5, 15),
                aspirate_rate=0.5,
                dispense_rate=0.5
            )
            
            # Slow tip withdrawal and drop
            slow_tip_withdrawal(p50m, firststrand_wells[col_idx][0], -2)
            p50m.drop_tip()
        
        # incubate 5 min with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=2.5*DRY_RUN)
        bead_mixing(p50m, 5, 30, speed=0.1)
        protocol.delay(minutes=2.5*DRY_RUN)

    if not skip_wash1: 
        # move plate to magnet for 2 min
        protocol.move_labware(labware=plate2_sample, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)
      
        #Remove and discard supernatant
        for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
            remove_sup(p50m, 20, 5, waste_well=-1)
            p50m.drop_tip()
    
        # two repeats in total
        for repeat in range(2):
            # move plate off magnet for washes
            for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
                if index < 2:
                    source = ethanol[0]
                elif index < 4:
                    source = ethanol[1]
                elif index < 6:
                    source = ethanol[2]


                pick_up_or_refill(p1000m)
                p1000m.aspirate(180, ethanol[0].bottom(clearance_reservoir))
                p1000m.dispense(180, column[0].bottom(1), rate = 0.5)
                p1000m.air_gap(20)  # prevent ethanol from dripping
                p1000m.drop_tip()

            # move plate to magnet for supernatant removal
            protocol.delay(minutes=1*BEAD_RUN)

            # remove supernatant completely
            for column in firststrand_wells[:NUM_COLUMNS]:
                remove_sup(p1000m,130,50,waste_well=repeat)
                
                if repeat == 1: 
                    p1000m.move_to(column[0].top())
                for clearance in [0.7, 0.4, 0.2, 0]:
                    loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0, z=0))
                    p1000m.aspirate(25, loc, rate= 0.2)
                p1000m.air_gap(20)
                p1000m.dispense(45, waste[1].top())
                p1000m.drop_tip()

        protocol.delay(seconds=30*BEAD_RUN)
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    if not skip_elution:
        for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
            pick_up_or_refill(p50m)
            p50m.aspirate(20, NFW.bottom(clearance_reservoir))
            
            # dispense Tris quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(20, loc, rate=2.5)
            wash_mixing(p50m, volume=16,reps=8)

            # side touch with blowout after last mix
            side_touch_w_blowout(p50m,column[0],pos=-5)       
            p50m.drop_tip()
        
        # incubate 10 min total with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=3*DRY_RUN)
        bead_mixing(p50m, 5, 15, speed=0.05)
        protocol.delay(minutes=3*DRY_RUN)
        bead_mixing(p50m, 5, 15, speed=0.05)
        protocol.delay(minutes=4*DRY_RUN)

        # transfer eluate to fresh wells
        for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
            pick_up_or_refill(p50m)
            p50m.move_to(column[0].top())
            p50m.air_gap(10)
            p50m.aspirate(20, column[0].bottom(4), rate=0.1)        
            protocol.delay(seconds=1)
            p50m.dispense(30, secondstrand_wells[index][0].bottom(1), rate=0.5)
            p50m.drop_tip()

    if not skip_2ndstrand:
        # distribute MM2 to required wells on the sample plate
        p50m.transfer(
            volume=30,
            source=MM3[0].bottom(1),
            dest=[col[0].bottom(1) for col in secondstrand_wells[:NUM_COLUMNS]],
            new_tip='always',
            aspirate_rate=0.5,
            dispense_rate=1,
            mix_after=(3, 10),  # Mix 3 times with 7.5µL after dispensing
            blow_out=True,
            blow_out_location='destination well',
            touch_tip=12    # Touch tip 12mm below the top of the well, check if this is appropriate
        )

        thermocycler.set_lid_temperature(100)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(94, hold_time_minutes=1*DRY_RUN, block_max_volume=50)
        thermocycler.set_block_temperature(50, hold_time_minutes=1*DRY_RUN, block_max_volume=50)
        thermocycler.set_block_temperature(65, hold_time_minutes=15*DRY_RUN, block_max_volume=50)
        thermocycler.set_block_temperature(20)
        thermocycler.open_lid()
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    if not skip_bead40:
        # Calculate mixing volumes for each bead column based on NUM_COLUMNS
        bead_mix_vols = []
        if NUM_COLUMNS <= 3:
            # Only use first column of beads
            bead_mix_vols.append(30 * NUM_COLUMNS)  # 50µl per sample column
            bead_mix_vols.append(0)  # Second column not used
        else:
            # Use both bead columns
            bead_mix_vols.append(90)  # First column always 150µl for 3 columns
            bead_mix_vols.append(30 * (NUM_COLUMNS - 3))  # Second column for remaining columns

        # Resuspend beads in each column with p1000
        for idx, (bead_col, mix_vol) in enumerate(zip(beads, bead_mix_vols)):
            if mix_vol > 0:  # Only process columns with beads
                p1000m.pick_up_tip()
                
                # Thorough mixing at different heights
                p1000m.mix(5, mix_vol * 0.8, bead_col[0].bottom(2))  # Bottom mixing
                p1000m.mix(5, mix_vol * 0.8, bead_col[0].bottom(5))  # Middle mixing
                
                # Final slow mixes
                for _ in range(3):
                    p1000m.aspirate(mix_vol, bead_col[0].bottom(2), rate=0.5)
                    p1000m.dispense(mix_vol, bead_col[0].bottom(5), rate=0.5)
                
                # Slow withdrawal
                slow_tip_withdrawal(p1000m, bead_col[0], -2)
                p1000m.drop_tip()

        # Distribute 17µl of beads to sample columns
        for col_idx in range(NUM_COLUMNS):
            bead_source = beads[0] if col_idx < 3 else beads[1]  # Switch source after 3 columns
            
            p50m.pick_up_tip()
            p50m.mix(3, 30, bead_source[0].bottom(2), rate=0.5)  # Gentle remix before aspiration
            
            # Transfer beads with mixing
            p50m.transfer(
                40,
                bead_source[0].bottom(2),
                secondstrand_wells[col_idx][0].bottom(2),
                new_tip='never',
                mix_after=(5, 50),
                aspirate_rate=0.5,
                dispense_rate=0.5
            )
            
            # Slow tip withdrawal and drop
            slow_tip_withdrawal(p50m, secondstrand_wells[col_idx][0], -2)
            p50m.drop_tip()
        
        # incubate 5 min with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=2.5*DRY_RUN)
        bead_mixing(p50m, 5, 50, speed=0.1)
        protocol.delay(minutes=2.5*DRY_RUN)

    if not skip_wash2: 
        # move plate to magnet for 2 min
        protocol.move_labware(labware=plate2_sample, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)
      
        #Remove and discard supernatant
        for index, column in enumerate(secondstrand_wells[:NUM_COLUMNS]):
            remove_sup(p50m, 20, 5, waste_well=-1)
            p50m.drop_tip()
    
        # two repeats in total
        for repeat in range(2):
            # move plate off magnet for washes
            for index, column in enumerate(secondstrand_wells[:NUM_COLUMNS]):
                if index < 2:
                    source = ethanol[0]
                elif index < 4:
                    source = ethanol[1]
                elif index < 6:
                    source = ethanol[2]


                pick_up_or_refill(p1000m)
                p1000m.aspirate(180, ethanol[0].bottom(clearance_reservoir))
                p1000m.dispense(180, column[0].bottom(1), rate = 0.5)
                p1000m.air_gap(20)  # prevent ethanol from dripping
                p1000m.drop_tip()

            # move plate to magnet for supernatant removal
            protocol.delay(minutes=1*BEAD_RUN)

            # remove supernatant completely
            for column in secondstrand_wells[:NUM_COLUMNS]:
                remove_sup(p1000m,130,50,waste_well=repeat)
                
                if repeat == 1: 
                    p1000m.move_to(column[0].top())
                for clearance in [0.7, 0.4, 0.2, 0]:
                    loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0, z=0))
                    p1000m.aspirate(25, loc, rate= 0.2)
                p1000m.air_gap(20)
                p1000m.dispense(45, waste[1].top())
                p1000m.drop_tip()

        protocol.delay(seconds=30*BEAD_RUN)
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    if not skip_finalelution:
        # TODO remove this section once tiprack exchange logic is finished or obsolete
        # move polyA plate to staging area and move eluate plate to deck
        # protocol.move_labware(labware=plate1_polyA, new_location="B4", use_gripper=True)
        # protocol.move_labware(labware=plate3_eluate, new_location="B2", use_gripper=True)

        # elute cDNA from beads
        for index, column in enumerate(secondstrand_wells[:NUM_COLUMNS]):
            pick_up_or_refill(p50m)
            p50m.aspirate(ELUTION_VOL, NFW.bottom(clearance_reservoir))
            
            # dispense quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(ELUTION_VOL, loc, rate=2.5)
            wash_mixing(p50m, volume=16,reps=8)

            # side touch with blowout after last mix
            side_touch_w_blowout(p50m,column[0],pos=-5)       
            p50m.drop_tip()
        
        # incubate 10 min total with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=3*DRY_RUN)
        bead_mixing(p50m, 5, 17, speed=0.05)
        protocol.delay(minutes=3*DRY_RUN)
        bead_mixing(p50m, 5, 17, speed=0.05)
        protocol.delay(minutes=4*DRY_RUN)

        # transfer eluate to fresh wells
        for index, column in enumerate(firststrand_wells[:NUM_COLUMNS]):
            pick_up_or_refill(p50m)
            p50m.move_to(column[0].top())
            p50m.air_gap(10)
            p50m.aspirate(20, column[0].bottom(4), rate=0.1)        
            protocol.delay(seconds=1)
            p50m.dispense(30, eluate_wells[index][0].bottom(1), rate=0.5)
            p50m.drop_tip()
