from opentrons import protocol_api
from opentrons import types
from opentrons.protocol_api.labware import OutOfTipsError

metadata = {
    "protocolName": "Reverse transcription and strand-switching (SQK-LSK114)",
    "description": "Based on the Nanopore direct cDNA sequencing protocol (DCS_9187_v114_revJ_12Dec2024)",
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
    DRY_RUN = 0.01      # use 0.01 to shorten wait times if it is dry run, otherwise 1
    BEAD_RUN = 0.01        # use 0.01 for testing without beads, otherwise 1 

    skip_firststrand = False     # Toggle when testing certain blocks, same below
    skip_65_5min = True
    skip_strandswitch = False
    skip_42_90min = True
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
    temp_block = protocol.load_module('temperature module gen2', 'D1')  
    temp_adapter = temp_block.load_adapter('opentrons_96_well_aluminum_block')     
    cold_plate = temp_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")

    ## COLUMN 2 
    mag_block = protocol.load_module('magneticBlockV1', 'A2')
    plate1_polyA = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "B2")
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    plate2_sample = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "D2")
    
    ## COLUMN 3 
    tiprack_50ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'A3')
    tiprack_50ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B3')
    tiprack_200ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C3')
    waste_chute = protocol.load_waste_chute()
    
    ## COLUMN 4 - staging area
    # tiprack_50ul_3 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'A4')
    # tiprack_50ul_4 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B4')
    # tiprack_200ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C4')
    # D4 is an empty slot to move tip boxes to
 
    #======== PIPETTES ========
    # load active tips for each pipette
    p50m = protocol.load_instrument('flex_8channel_50', 'right')
    p1000m = protocol.load_instrument('flex_8channel_1000', 'left')
    
    # define "active" and "staging" tip boxes for each pipette
    p50m.tip_racks = [tiprack_50ul_1, tiprack_50ul_2]
    # p50_staging = [tiprack_50ul_3, tiprack_50ul_4]
    
    p1000m.tip_racks = [tiprack_200ul_1]
    # p1000_staging = [tiprack_200ul_2]

    #======== REAGENT WELLS ========
    # define well locations for where samples are moved from/to 
    polyA_wells = plate1_polyA.columns()[6 : 6 + NUM_COLUMNS]  # oligo dT-bead eluates from polyA protocol, old plate columns 7-12)
    firststrand_wells = plate2_sample.columns()[:NUM_COLUMNS]  # wells for the first strand synthesis
    secondstrand_wells = plate2_sample.columns()[6 : 6 + NUM_COLUMNS] # wells for the second strand synthesis

    # cold plate (for reagents and the eluate for this protocol)
    MM1 = cold_plate.columns_by_name()['1']
    MM2 = cold_plate.columns_by_name()['2']
    MM3 = cold_plate.columns_by_name()['3']
    beads = [cold_plate.columns_by_name()[name] for name in ['4', '5']]
    rnase = cold_plate.columns_by_name()['6']
    eluate_wells = cold_plate.columns()[6 : 6 + NUM_COLUMNS] # wells for the final eluate

    # beads (2X 47 uL per sample)
    for index, column in enumerate(beads):
        if not index:
            num = NUM_COLUMNS if NUM_COLUMNS <= 3 else 3    # only 1 column for up to 3 samples
        else:
            num = NUM_COLUMNS - 3 if NUM_COLUMNS > 3 else 0 # 2 columns for more than 3 samples
        # # sets volume needed in each column
        # # note: I don't know where this is used if at all
        # if num:
        #     column[0].liq_vol = 50*num + deadvol_plate
        # else:
        #     column[0].liq_vol = 0

    # reagent reservoir
    ethanol = [reservoir.wells_by_name()[well] for well in ['A1', 'A2', 'A3']]
    # for well in ethanol: well.liq_vol = 4*200*2*p1000m.channels + deadvol_reservoir
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

    #======== RUN SETUP ========
    # HELPER FUNCTIONS 
    # FUNCTION 1: PIPETTE PICK UP
    def pick_up_or_refill(pipette):
        """Pick up a tip or pause for replacement if needed."""
        try:
            pipette.pick_up_tip()
        except OutOfTipsError:
            protocol.pause(
             """Please Refill the {} Tip Boxes
                and Empty the Tip Waste.""".format(pipette.mount))
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
    def side_touch_w_blowout(pipette,well,pos=-1,xoffset=0):
            pipette.default_speed /= 40
            pipette.move_to(well.top(pos).move(types.Point(x=(well.diameter / 2)-xoffset, y=0, z=0))) 
            pipette.blow_out()
            protocol.delay(seconds=1.5)
            pipette.default_speed *= 2
            pipette.move_to(well.top())
            pipette.default_speed *= 20

    # FUNCTION 4: MIX WITH SLOW WITHDRAWAL
    def slow_mixing(pipette, columns, reps, vol, speed=0.5):
        for index, column in enumerate(columns):
            if not pipette.has_tip: 
                pick_up_or_refill(pipette)
            
            # slow mixes
            pipette.mix(reps, vol, column[0].bottom(2), rate=speed)
            slow_tip_withdrawal(pipette, column[0], -2)
            
            # side touch with blowout after last mix
            side_touch_w_blowout(pipette,column[0],pos=-5)                  
            pipette.drop_tip()

    # FUNCTION 6: PELLET MIXING (for when beads are in a pellet)
    def pellet_mixing(pipette, column, volume,reps):
        
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
        pipette.mix(reps, volume, column[0].bottom(2), rate=0.5)
        slow_tip_withdrawal(pipette, column[0], -2)
        
    # FUNCTION 5: REMOVE SUPERNATANT
    def remove_sup(pipette,column,volume1,volume2,waste_well):       
        pick_up_or_refill(pipette)
        
        pipette.move_to(column[0].top())
        pipette.air_gap(20)
        pipette.aspirate(volume1, column[0].bottom(4), rate=0.2)        
        pipette.aspirate(volume2, column[0].bottom(clearance_bead_pellet), rate=0.1)
        protocol.delay(seconds=1)

        pipette.dispense(volume1+volume2+20, waste[waste_well].top(), rate=2)
        protocol.delay(seconds=1)
        pipette.blow_out()

    def add_beads(destination_wells, bead_volume, mix_volume_multiplier, final_mix_volume):
        """
        Perform bead mixing and distribution steps
        
        Args:
            NUM_COLUMNS: Number of columns to process
            beads: Bead source wells
            destination_wells: Destination wells for bead transfer
            bead_volume: Volume of beads to transfer (17 or 40 µL)
            mix_volume_multiplier: Volume multiplier for bead mixing (50 or 30)
            final_mix_volume: Final mixing volume after transfer (15 or 45)
        """
        # Calculate mixing volumes for each bead column
        bead_mix_vols = []
        if NUM_COLUMNS <= 3:
            bead_mix_vols.append(mix_volume_multiplier * NUM_COLUMNS)
            bead_mix_vols.append(0)  # Second column not used
        else:
            bead_mix_vols.append(mix_volume_multiplier * 3)  # First column for 3 samples
            bead_mix_vols.append(mix_volume_multiplier * (NUM_COLUMNS - 3))  # Second column for remaining

        # Resuspend beads in source columns
        for bead_col, mix_vol in zip(beads, bead_mix_vols):
            if mix_vol > 0:
                pick_up_or_refill(p1000m)
                pellet_mixing(pipette=p1000m, column=bead_col, volume=mix_vol,reps=5)
                side_touch_w_blowout(p1000m,bead_col[0],pos=-2)
                p1000m.drop_tip()

        # Distribute beads to destination wells
        for index, column in enumerate(destination_wells):
            bead_source = beads[0] if index < 3 else beads[1]
            
            # quick mix before transfer
            pick_up_or_refill(p50m)
            p50m.mix(3, mix_volume_multiplier, bead_source[0].bottom(2), rate=1)
        
            # transfer beads and mix after
            p50m.aspirate(bead_volume, bead_source[0].bottom(1), rate=0.5)
            p50m.dispense(bead_volume, destination_wells[index][0].bottom(1), rate=1)
            protocol.delay(seconds=0.5)
            p50m.mix(5, final_mix_volume, destination_wells[index][0].bottom(1), 1)
            slow_tip_withdrawal(p50m, destination_wells[index][0], -2)
            side_touch_w_blowout(p50m,destination_wells[index][0],pos=-8,xoffset=1)
            p50m.drop_tip()
    
    # FUNCTION 7: ETHANOL RINSE (washing on the magnet block)
    def ethanol_rinse(sup_pip, sup_vol, plate, wells, waste):
        # move plate to magnet for 2 min
        protocol.move_labware(labware=plate, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)
    
        #Remove and discard supernatant
        for index, column in enumerate(wells):
            remove_sup(sup_pip, column, sup_vol-5, 5, waste_well=-1)
            sup_pip.drop_tip()

        # two repeats in total
        for repeat in range(2):
            # move plate off magnet for washes
            for index, column in enumerate(wells):
                if index < 2:
                    source = ethanol[0]
                elif index < 4:
                    source = ethanol[1]
                elif index < 6:
                    source = ethanol[2]

                pick_up_or_refill(p1000m)
                p1000m.aspirate(180, ethanol[0].bottom(clearance_reservoir))
                p1000m.dispense(180, column[0].bottom(5), rate = 0.1)
                p1000m.air_gap(20)  # prevent ethanol from dripping
                p1000m.drop_tip()

            # move plate to magnet for supernatant removal
            protocol.delay(minutes=1*BEAD_RUN)

            # remove supernatant completely
            for index, column in enumerate(wells):
                remove_sup(p1000m,column,130,50,waste_well=repeat)
                # try to not do drop tip
                if repeat == 1: 
                    p1000m.move_to(column[0].top())
                for clearance in [0.7, 0.4, 0.2, 0]:
                    loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0, z=0))
                    p1000m.aspirate(25, loc, rate= 0.2)
                p1000m.air_gap(20)
                p1000m.dispense(45, waste[1].top())
                p1000m.blow_out()
                p1000m.drop_tip()

        protocol.delay(seconds=30*BEAD_RUN)
        protocol.move_labware(labware=plate, new_location="D2", use_gripper=True)
        
    #======== STEP BY STEP PROTOCOL ========
    # Set the temperature of the cooler block and thermo cycler
    temp_block.set_temperature(celsius=4)  
    thermocycler.set_block_temperature(temperature=20)
    
    ## 4. Reverse transcription and strand-switching
    # STEP 4.3 - mix RNA input and VN primer mix
    if not skip_firststrand:
        
        # distribute 3.5µl MM1 to required wells on the sample plate
        transfer_vol = 3.5
        dead_vol = 2
        destinations = [col[0].bottom(1) for col in firststrand_wells]

        pick_up_or_refill(p50m)
        p50m.aspirate((transfer_vol * len(destinations)) + dead_vol, 
                      MM1[0].bottom(1), rate=0.5)
        for d in destinations:
            p50m.dispense(transfer_vol, d, rate=0.3)
            protocol.delay(seconds=0.5)
        p50m.drop_tip()
        
        # Transfer polyA samples to sample plate and mix
        polyA_vol = 7.5
        for src_col, dest_col in zip(polyA_wells, firststrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(polyA_vol, src_col[0].bottom(1), rate=0.5)
            p50m.dispense(polyA_vol, dest_col[0].bottom(1), rate=1.0)
            p50m.mix(3, polyA_vol, dest_col[0].bottom(1), 0.5)
            side_touch_w_blowout(p50m,dest_col[0],pos=-11,xoffset=2)
            p50m.drop_tip()

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
        for index, column in enumerate(firststrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(9, MM2[0].bottom(1), rate=0.5)
            p50m.dispense(9, column[0].bottom(1), rate=0.5)
            p50m.mix(3, 15, column[0].bottom(1), rate=0.5)
            side_touch_w_blowout(p50m, column[0], pos=-9, xoffset=1)
            p50m.drop_tip()

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
        RNAse_vol=1
        
        for index, column in enumerate(firststrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(RNAse_vol, rnase[0].bottom(1), rate=0.5)
            p50m.dispense(RNAse_vol, column[0].bottom(2.5), rate=0.5)
            protocol.delay(seconds=0.5)
            p50m.mix(3, 17, column[0].bottom(1), 0.5)
            side_touch_w_blowout(p50m,column[0],pos=-9,xoffset=1)
            p50m.drop_tip()

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
        add_beads(
            destination_wells=firststrand_wells,
            bead_volume=17,
            mix_volume_multiplier=50,
            final_mix_volume=30)
        
        # incubate 5 min with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=2.5*DRY_RUN)
        slow_mixing(p50m, firststrand_wells, reps=5, vol=30, speed=0.5)
        protocol.delay(minutes=2.5*DRY_RUN)

    if not skip_wash1: 
        ethanol_rinse(sup_pip=p50m, sup_vol=30, plate=plate2_sample, wells=firststrand_wells, waste=waste)

    if not skip_elution:
        for index, column in enumerate(firststrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(20, NFW.bottom(clearance_reservoir))
            
            # dispense Tris quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(20, loc, rate=2.5)
            pellet_mixing(p50m, column, volume=16,reps=8)

            # side touch with blowout after last mix
            side_touch_w_blowout(p50m,column[0],pos=-5,xoffset=1)       
            p50m.drop_tip()
        
        # incubate 10 min total with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=3*DRY_RUN)
        slow_mixing(p50m, firststrand_wells, reps=5, vol=15, speed=0.3)
        protocol.delay(minutes=3*DRY_RUN)
        slow_mixing(p50m, firststrand_wells, reps=5, vol=15, speed=0.3)
        protocol.delay(minutes=4*DRY_RUN)

        # place on magnet to pellet beads
        protocol.move_labware(labware=plate2_sample, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=1.5*BEAD_RUN)

        # transfer eluate to fresh wells
        for index, column in enumerate(firststrand_wells):
            pick_up_or_refill(p50m)
            p50m.move_to(column[0].top())
            p50m.air_gap(10)
            p50m.aspirate(20, column[0].bottom(clearance_bead_pellet), rate=0.1)        
            protocol.delay(seconds=1)
            p50m.dispense(30, secondstrand_wells[index][0].bottom(1), rate=0.2)
            side_touch_w_blowout(p50m,secondstrand_wells[index][0],pos=-8,xoffset=1)
            p50m.drop_tip()
        
        # move plate off magnet
        protocol.move_labware(labware=plate2_sample, new_location="D2", use_gripper=True)

    if not skip_2ndstrand:
        # distribute MM3 to required wells on the sample plate
        for index, column in enumerate(secondstrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(30, MM3[0].bottom(1), rate=1)
            p50m.dispense(30, column[0].bottom(1), rate=1)
            p50m.mix(3, 40, column[0].bottom(1), rate=0.5)
            side_touch_w_blowout(p50m, column[0], pos=-7, xoffset=1)  # Check if -9 offset is correct
            p50m.drop_tip()
        
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
        add_beads(
            destination_wells=secondstrand_wells,
            bead_volume=40,
            mix_volume_multiplier=30,
            final_mix_volume=45)
        
        # incubate 5 min with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=2.5*DRY_RUN)
        slow_mixing(p50m, secondstrand_wells, reps=5, vol=50, speed=0.5)
        protocol.delay(minutes=2.5*DRY_RUN)

    if not skip_wash2: 
        ethanol_rinse(sup_pip=p1000m, sup_vol=90, plate=plate2_sample, wells=secondstrand_wells, waste=waste)
    # TODO change back to 80 after figuring out the volume tracking error

    if not skip_finalelution:
        # TODO remove this section once tiprack exchange logic is finished or obsolete
        # move polyA plate to staging area and move eluate plate to deck
        # protocol.move_labware(labware=plate1_polyA, new_location="B4", use_gripper=True)

        # elute cDNA from beads
        for index, column in enumerate(secondstrand_wells):
            pick_up_or_refill(p50m)
            p50m.aspirate(ELUTION_VOL, NFW.bottom(clearance_reservoir))
            
            # dispense quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(ELUTION_VOL, loc, rate=2.5)
            pellet_mixing(p50m, column, volume=16,reps=8)

            # side touch with blowout after last mix
            side_touch_w_blowout(p50m,column[0],pos=-5)       
            p50m.drop_tip()
        
        # incubate 10 min total with slow mixing (Hula mixer replacement)  
        protocol.delay(minutes=3*DRY_RUN)
        slow_mixing(p50m, secondstrand_wells, reps=5, vol=17, speed=0.3)
        protocol.delay(minutes=3*DRY_RUN)
        slow_mixing(p50m, secondstrand_wells, reps=5, vol=17, speed=0.3)
        protocol.delay(minutes=4*DRY_RUN)

        protocol.move_labware(labware=plate2_sample, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)

        # transfer eluate to fresh wells
        for index, column in enumerate(secondstrand_wells):
            pick_up_or_refill(p50m)
            p50m.move_to(column[0].top())
            p50m.air_gap(10)
            p50m.aspirate(20, column[0].bottom(4), rate=0.2)        
            protocol.delay(seconds=1)
            p50m.dispense(30, eluate_wells[index][0].bottom(1), rate=0.5)
            p50m.drop_tip()
