from opentrons import protocol_api
from opentrons import types
from opentrons.protocol_api.labware import OutOfTipsError

metadata = {
    "protocolName": "NEBNExt High Input Poly(A) mRNA isolation",
    "description": "Protocol to work with the NEBNext poly(A) mRNA isolation kit NEB #E3370S",
    "author": "Haiyin Liu"
    }

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.18",
}

def run(protocol: protocol_api.ProtocolContext):
    #======== PARAMETERS ========
    # SAMPLE PARAMETERS
    NUM_SAMPLES = 48  # Define the number of samples (up to 48)
    NUM_COLUMNS = (NUM_SAMPLES + 7) // 8  # Calculate number of columns for 96-well plates
    
    ELUTION_VOL = 15    # µl of Tris to resuspend the beads in at the last elution
    COLLECTION_VOL = 12    # µl of Tris to collect from the beads at the final magnet step

    # TESTING PARAMETERS
    DRY_RUN = 1      # use 0.01 to shorten wait times if it is dry run, otherwise 1
    BEAD_RUN = 1        # use 0.01 for testing without beads

    skip_mixbeadsandrna = False     # Toggle when testing certain blocks, same below
    skip_thermocycler1 = False
    skip_beadmix = False
    skip_supremoval = False
    skip_washes = False
    skip_1stelution = False
    skip_thermocycler2 = False
    skip_rebinding = False
    skip_supremoval2 = False
    skip_finalwash = False
    skip_washremoval = False
    skip_2ndelution = False
    skip_heat_collect = False

    # VOLUME AND DISTANCE SETTINGS
    deadvol_reservoir = 1500
    deadvol_plate = 10
    clearance_reservoir = 2
    clearance_bead_pellet = 1.5
    clearance_beadresuspension = 3
    offset_x = 1

    #======== DECK SETUP ========
    # COLUMN 1 
    thermocycler = protocol.load_module('thermocycler module gen2')
    temp_block = protocol.load_module('temperature module gen2', 'D1')
    temp_adapter = temp_block.load_adapter('opentrons_96_well_aluminum_block')     
    cold_plate = temp_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    
    # COLUMN 2 
    mag_block = protocol.load_module('magneticBlockV1', 'A2')
    tiprack_50ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B2')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    sample_plate = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "D2")
     
    # COLUMN 3 
    tiprack_200ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'A3')
    tiprack_200ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'B3')
    tiprack_200ul_3 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C3')
    waste_chute = protocol.load_waste_chute()
     
    # COLUMN 4 - staging area
    tiprack_200ul_4 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'A4')
    tiprack_200ul_5 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'B4')
    tiprack_200ul_6 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C4')
    # D4 is an empty slot to move tip boxes to

    #======== PIPETTES ========
    p50m = protocol.load_instrument('flex_8channel_50', 'right')
    p1000m = protocol.load_instrument('flex_8channel_1000', 'left')

    p50m.tip_racks = [tiprack_50ul_1]

    # p1000m pipette has "active" and "staging" tip boxes for swapping
    p1000m.tip_racks = [tiprack_200ul_1,tiprack_200ul_2,tiprack_200ul_3]
    p1000_staging = [tiprack_200ul_4,tiprack_200ul_5,tiprack_200ul_6]

    # tracking if the staging tip boxes have been swapped in
    has_swapped = {
        "p1000_multi_flex": False
    } 

    #======== REAGENT WELLS ========
    # sample plate 
    sample_wells = sample_plate.columns()[:NUM_COLUMNS]
    
    # cold plate (for reagents and the eluate)
    rna_binding_buffer = cold_plate.columns()[: (1 if NUM_COLUMNS <= 3 else 2)] # column #1 (+ #2 for more than 3 samples)
    tris = cold_plate.columns()[3 : (4 if NUM_COLUMNS <= 3 else 5)]             # column #4 (+ #5 for more than 3 samples)
    finaleluate_wells = cold_plate.columns()[6 : 6 + NUM_COLUMNS]               # columns #6-12 for eluate

    # binding buffer (2X 50 uL per sample)
    for index, column in enumerate(rna_binding_buffer):
        if not index:
            num = NUM_COLUMNS if NUM_COLUMNS <= 3 else 3    # only 1 column for up to 3 samples
        else:
            num = NUM_COLUMNS - 3 if NUM_COLUMNS > 3 else 0 # 2 columns for more than 3 samples

    # reagent reservoir
    ### TODO check if volume correct - wash buffer 180 uL per sample - for three washes
    wash_buffer = [reservoir.wells_by_name()[well] for well in ['A1','A2','A3']]
      
    waste = [reservoir.wells_by_name()[well] for well in ['A10', 'A11', 'A12']]

    #======== DEFINING LIQUIDS ========   
    sample_liquid = protocol.define_liquid(name="Samples",
                                          description="50 µL sample + 50 µL 2X binding buffer",
                                          display_color="#C0C0C0") # Silver
    sample_vol = 100    # 50 µL RNA sample with 50 µL 2x binding buffer added just prior

    wash_buffer_liquid = protocol.define_liquid(name="wash buffer", 
                                               description="wash buffer, 180µl * sample * number of washes", 
                                               display_color="#0000FF") # Blue
    wash_buffer_vol = 180 * NUM_SAMPLES * 3 + deadvol_reservoir

    tris_liquid = protocol.define_liquid(name="Tris buffer",
                                        description="50 µL per sample * sample count",
                                        display_color="#FFA500") # Orange
    tris_vol = 50 * NUM_COLUMNS/2 + deadvol_plate

    binding_buffer_liquid = protocol.define_liquid(name="2X RNA binding buffer",
                                        description="50 µL per sample * sample count",
                                        display_color="#008000") # Green
    binding_buffer_vol = 65 * NUM_COLUMNS/2 + deadvol_plate

    eluate_liquid = protocol.define_liquid(name="Eluate",
                                           description="Eluted RNA after cleanup",
                                           display_color="#0000FF") # Blue
    eluate_vol = 15

    #======== LOADING LIQUIDS ========
    # sample plate
    for well in sample_plate.wells()[:NUM_SAMPLES]:
        well.load_liquid(liquid=sample_liquid, volume=sample_vol)

    # reservoir
    for column in reservoir.columns()[0:2]:
        for well in column: 
            well.load_liquid(liquid=wash_buffer_liquid, volume=wash_buffer_vol)

    # cold plate
    # dictionary for mapping liquids on the plate 
    liquid_mapping = {
        # Binding buffer goes in columns "1" and "2"
        '1': (binding_buffer_liquid, binding_buffer_vol),
        '2': (binding_buffer_liquid, binding_buffer_vol),
        
        # Tris buffer goes in columns "4" and "5"
        '4': (tris_liquid, tris_vol),
        '5': (tris_liquid, tris_vol),
        
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
    def pick_up_or_swap(pipette):
        """
        Try picking up a tip. If out of tips:
        - p50m: Pause for manual refill.
        - p1000m: If not swapped yet, swap 3 active racks, otherwise do full manual refill.
        """
        pip_name = pipette.name  # e.g. 'p1000_multi_flex' or 'p50_multi_flex'
        try:
            pipette.pick_up_tip()
        except OutOfTipsError:
            protocol.comment(f"{pip_name} is out of tips.")
            
            # --- p50m logic: always manual refill ---
            if pipette == p50m:
                protocol.pause("The 50µL tip box is empty. Please replace it manually.")
                pipette.reset_tipracks()
                pipette.pick_up_tip()
                return

            # --- p1000m logic: triple swap for the 3 active racks, if not done yet ---
            if has_swapped[pip_name]:
                # Already did our triple-swap => must do a full manual replacement now
                protocol.pause("No more 200µL tip boxes left. Please replace racks manually.")
                pipette.reset_tipracks()
                pipette.pick_up_tip()
                return
            else:
                # Not swapped yet => do the 3 triple-moves at once
                protocol.comment(f"{pip_name} is out of tips. Swapping racks from staging...")

                # SWAP #1:
                # tiprack_200ul_1 -> D4, tiprack_200ul_4 -> A3, tiprack_200ul_1 -> A4
                protocol.move_labware(tiprack_200ul_1, "D4", use_gripper=True)
                protocol.move_labware(tiprack_200ul_4, "A3", use_gripper=True)
                protocol.move_labware(tiprack_200ul_1, "A4", use_gripper=True)
                if tiprack_200ul_1 in pipette.tip_racks:
                    idx = pipette.tip_racks.index(tiprack_200ul_1)
                    pipette.tip_racks[idx] = tiprack_200ul_4

                # SWAP #2:
                # tiprack_200ul_2 -> D4, tiprack_200ul_5 -> B3, tiprack_200ul_2 -> B4
                protocol.move_labware(tiprack_200ul_2, "D4", use_gripper=True)
                protocol.move_labware(tiprack_200ul_5, "B3", use_gripper=True)
                protocol.move_labware(tiprack_200ul_2, "B4", use_gripper=True)
                if tiprack_200ul_2 in pipette.tip_racks:
                    idx = pipette.tip_racks.index(tiprack_200ul_2)
                    pipette.tip_racks[idx] = tiprack_200ul_5

                # SWAP #3:
                # tiprack_200ul_3 -> D4, tiprack_200ul_6 -> C3, tiprack_200ul_3 -> C4
                protocol.move_labware(tiprack_200ul_3, "D4", use_gripper=True)
                protocol.move_labware(tiprack_200ul_6, "C3", use_gripper=True)
                protocol.move_labware(tiprack_200ul_3, "C4", use_gripper=True)
                if tiprack_200ul_3 in pipette.tip_racks:
                    idx = pipette.tip_racks.index(tiprack_200ul_3)
                    pipette.tip_racks[idx] = tiprack_200ul_6

                # Mark that we've done our triple-swap
                has_swapped[pip_name] = True
                
                # Refresh tip tracking & pick up
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
    def bead_mixing(reps, vol, speed=0.5):
        for index, column in enumerate(sample_wells):
            if not p1000m.has_tip: 
                pick_up_or_swap(p1000m)
            
            # slow mixes
            p1000m.mix(reps, vol, column[0].bottom(3), rate=speed)
            slow_tip_withdrawal(p1000m, column[0], -2)
            
            # side touch with blowout after last mix
            side_touch_w_blowout(p1000m,column[0],pos=-1)                  
            p1000m.drop_tip()

    # FUNCTION 5: REMOVE SUPERNATANT
    def remove_sup(volume1,volume2,waste_well):       
        pick_up_or_swap(p1000m)
        
        p1000m.move_to(column[0].top())
        p1000m.air_gap(20)
        p1000m.aspirate(volume1, column[0].bottom(4), rate=0.1)        
        p1000m.aspirate(volume2, column[0].bottom(clearance_bead_pellet), rate=0.02)
        protocol.delay(seconds=1)

        p1000m.dispense(volume1+volume2+20, waste[waste_well].top(), rate=2)
        protocol.delay(seconds=1)
        p1000m.blow_out()

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
        

    # FUNCTION X: LIQUID HEIGHT IN A WELL
    # legacy code, I have not included this as we don't seem to need to calculate the height
    # TODO figure out if this is needed, and if so, import math
    #  def liq_height(well, effective_diameter=None):
    #     if well.diameter:
    #         if effective_diameter:
    #             radius = effective_diameter / 2
    #         else:
    #             radius = well.diameter / 2
    #         csa = math.pi*(radius**2)
    #     else:
    #         csa = well.length*well.width
    #     return well.liq_vol / csa

    # STEP-BY-STEP PROTOCOL
    # Set the temperature of the cooler block and thermo cycler
    temp_block.set_temperature(celsius=4)  
    thermocycler.set_block_temperature(temperature=20)
    
    # STEP 11: Mix beads and RNA sample 10 times
    protocol.comment("Step 11: Mixing beads/RNA sample by pipetting.")
    
    if not skip_mixbeadsandrna:
        
        bead_mixing(10, 80) # 10 mixes with 80µl 

    # STEP 12: Heat for RNA denaturation and binding
    protocol.comment("Step 12: Heating samples in thermocycler for RNA denaturation and binding.")
    
    if not skip_thermocycler1:
        thermocycler.set_lid_temperature(75)
        thermocycler.open_lid()
        protocol.move_labware(labware=sample_plate, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(65, hold_time_minutes=5*DRY_RUN, block_max_volume=100)
        thermocycler.set_block_temperature(4, hold_time_minutes=1*DRY_RUN)
        thermocycler.open_lid()
        
        # STEP 13: Remove tubes from thermocycler
        protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)
        thermocycler.set_block_temperature(temperature=20)

    # STEP 14-17: Resuspend beads on the bench, allow the RNA to bind to the beads, repeat 1x
    protocol.comment("Step 14-17: Pipette up and down slowly ten times to mix, bind for 5 min, repeat 1x")
    
    if not skip_beadmix:

        for repeat in range(2):
            
            bead_mixing(10, 80) # 10 mixes with 80µl 
            protocol.delay(minutes=5*DRY_RUN)

    # STEP 18: move plate to magnet for 2 min to separate the poly(A) RNA
        protocol.comment("Step 18: Separate bead-bound RNA from solution.")

        protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)

    # STEP 19: Remove and discard supernatant
    protocol.comment("Step 19: Remove and discard supernatant.")
    if not skip_supremoval:
        for index, column in enumerate(sample_wells):
            remove_sup(75,75,waste_well=-1)
            p1000m.air_gap(20)
            p1000m.drop_tip()
    
    # STEP 20-25: Add wash buffer and wash twice
    protocol.comment("Step 20-25: Take off magnet, wash beads and mix 10 times. Repeat")

    if not skip_washes:

        # washing twice total
        for repeat in range(2):
            
            # move plate off magnet for washes
            protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

            # dispense wash buffer to all bead pellets first to avoid beads drying
            pick_up_or_swap(p1000m)
            source = wash_buffer[0]

            for index, column in enumerate(sample_wells):
                if index < 2:
                    source = wash_buffer[0]
                elif index < 4:
                    source = wash_buffer[1]
                elif index < 6:
                    source = wash_buffer[2]

                p1000m.aspirate(180, wash_buffer[0].bottom(clearance_reservoir)) # ht instead of 2 if using
                p1000m.dispense(180, column[0].top(), rate = 0.8)
                # protocol.delay(seconds=1)
                # p1000m.blow_out()

            # mix all
            for index, column in enumerate(sample_wells):
                if not p1000m.has_tip:
                    pick_up_or_swap(p1000m)
                wash_mixing(p1000m, volume=144,reps=8)
                p1000m.drop_tip()

            # move plate to magnet for supernatant removal
            protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
            protocol.delay(minutes=2*BEAD_RUN)

            # remove supernatant
            for column in sample_wells:
                remove_sup(130,50,waste_well=repeat)
                p1000m.air_gap(20)
                p1000m.drop_tip()
          
        protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 26: adding Tris for elution
    protocol.comment("Step 26: First elution with Tris, mix 10 times.")

    if not skip_1stelution:
        source = tris[0][0]

        for index, column in enumerate(sample_wells):
            if index > 2:
                source = tris[1][0]

            pick_up_or_swap(p50m)
            p50m.aspirate(50, source.bottom(1))
            
            # dispense Tris quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(50, loc, rate=2.5)
            wash_mixing(p50m, volume=43,reps=8)

            # side touch with blowout after last mix
            side_touch_w_blowout(p50m,column[0],pos=-1)       
            p50m.drop_tip()

    # STEP 27: First elution in thermocycler
    protocol.comment("Step 27: Heating samples in thermocycler for first elution from beads.")
    
    if not skip_thermocycler2:
        thermocycler.set_lid_temperature(90)
        thermocycler.open_lid()
        protocol.move_labware(labware=sample_plate, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(80, hold_time_minutes=2*DRY_RUN, block_max_volume=50)
        thermocycler.set_block_temperature(25)
        thermocycler.open_lid()
        # STEP 28: Remove tubes from thermocycler 
        protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)
    
    # STEP 29: Add 50 μl of RNA Binding Buffer (2X), mix and re-bind, repeat 1x
    
    protocol.comment("Step 29-32: Add RNA binding buffer and mix, re-bind, repeat 1x")
    
    if not skip_rebinding:
        source = rna_binding_buffer[0][0]

        for index, column in enumerate(sample_wells):
            # use second column of Binding buffer for sample columns >3
            if index > 2:
                source = rna_binding_buffer[1][0]
            
            # add 2x Binding buffer
            pick_up_or_swap(p1000m)
            p1000m.aspirate(50, source.bottom(1))
            p1000m.dispense(50, column[0].bottom(2))
            
            # slow mixing
            bead_mixing(10, 80, speed=0.3) # 10 mixes with 80µl
        
        # rebinding #1
        protocol.delay(minutes=5*DRY_RUN)
        
        # slow mixing and rebinding #2
        for column in sample_wells:
            bead_mixing(10, 80, speed=0.3) # 10 mixes with 80µl
        
        protocol.delay(minutes=5*DRY_RUN)
    
    # STEP 33: move plate to magnet for 2 min to separate the poly(A) RNA
        protocol.comment("Step 33: Second binding - separate RNA from solution.")
        
        protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)

        # STEP 34-35: Remove and discard supernatant, move off magnet
        protocol.comment("Step 34: Remove and discard supernatant.")
        if not skip_supremoval2:
            for column in sample_wells:
                remove_sup(75,75,waste_well=-1)
                p1000m.air_gap(20)
                p1000m.drop_tip()

        protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 36: Add final wash buffer and mix
    protocol.comment("Step 36: Wash the beads once with 200 μl of Wash Buffer, mix thoroughly.")

    if not skip_finalwash:
        for index, column in enumerate(sample_wells):

            pick_up_or_swap(p1000m)

            ### TODO this section tracks the liquid height in the reservoir, check if needd 
            # wash_buffer[-1].liq_vol -= 180*p1000m.channels
            # ht = liq_height(
            #   wash_buffer[-1]) - 3 if liq_height(
            #   wash_buffer[-1]) - 3 > 1 else 1
            p1000m.aspirate(180, wash_buffer[0].bottom(2))
            p1000m.dispense(180, column[0].bottom(2)) # ht instead of 2 if using
            wash_mixing(p1000m, volume=144,reps=8)
            p1000m.drop_tip()

        # STEP 37: Place tubes on magnetic rack
        protocol.comment("STEP 37: Place the tubes on the magnetic rack for 2 minutes")
        protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*BEAD_RUN)

    # STEP 38-39: Remove wash buffer completely and move plate off magnet
    protocol.comment("STEP 38-39: Remove and discard wash buffer, remove from magnet")

    # TODO Double check the top(), bottom(2) and rate parameters
    if not skip_washremoval:
        for column in sample_wells:
            pick_up_or_swap(p1000m)
            
            p1000m.move_to(column[0].top())
            p1000m.air_gap(20)
            p1000m.aspirate(110, column[0].bottom(4), rate=0.2)
            
            p1000m.aspirate(
                70, column[0].bottom(1).move(types.Point(x=offset_x, y=0, z=0)), rate=0.2)
            
            p1000m.dispense(200, waste[1].top(), rate=2)
            protocol.delay(seconds=1)
            p1000m.blow_out()

            # loop to move closer to the bottom
            p1000m.move_to(column[0].top())
            for clearance in [0.7, 0.4, 0.2, 0]:
                loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0, z=0))
                p1000m.aspirate(25, loc, rate= 0.2)
            p1000m.air_gap(20)
            p1000m.dispense(45, waste[1].top(), rate=2)
            p1000m.drop_tip()
      
        protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 40: Last elution in thermocycler
    if not skip_2ndelution:
        protocol.comment("STEP 40: Final elution with Tris, mix and heat in thermocycler.")

        source = tris[0][0]
        
        for index, column in enumerate(sample_wells):
            if index > 2:
                source = tris[1][0]
        
            pick_up_or_swap(p50m)
            p50m.aspirate(ELUTION_VOL, source.bottom(1))
            
            # dispense Tris quickly on the side of pellet and mix
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p50m.dispense(ELUTION_VOL, loc, rate=1)
            
            # side mixes to resuspend beads
            locations = [
                types.Point(x=offset_x, y=0, z=0),
                types.Point(x=-offset_x, y=0, z=0),
                types.Point(x=0, y=offset_x, z=0),
                types.Point(x=0, y=-offset_x, z=0)]

            for point in locations:
                loc = column[0].bottom(1.5).move(point)
                p50m.mix(3, 7, loc, rate=0.5)
            
            # slower mixes with larger volume
            p50m.mix(6, 12, column[0].bottom(3), rate=0.2)
            slow_tip_withdrawal(p50m, column[0], -2)
            p50m.drop_tip()
            
    if not skip_heat_collect:
        thermocycler.set_lid_temperature(90)
        thermocycler.open_lid()
        protocol.move_labware(labware=sample_plate, new_location=thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(80, hold_time_minutes=2*DRY_RUN, block_max_volume=50)
        thermocycler.set_block_temperature(25)
        thermocycler.open_lid()
        
        # STEP 41: place tubes on magnetic rack
        protocol.comment("STEP 41: Put on maget for 2 minutes")
        
        protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=2*DRY_RUN)
    
        # STEP 42: collect purified mRNA
        protocol.comment("STEP 42: Transfer supernatant into clean wells.")
        
        for index, source_column in enumerate(sample_wells):
            pick_up_or_swap(p50m)
            
            source_well = source_column[0]
            target_well = finaleluate_wells[index][0]
            
            p50m.move_to(source_well.top())
            p50m.air_gap(10)
            p50m.aspirate(COLLECTION_VOL, source_well.bottom(1), rate=0.2)
            protocol.delay(seconds=1)

            p50m.dispense(COLLECTION_VOL+5, target_well.bottom(1), rate=0.2)
            protocol.delay(seconds=2)
            p50m.blow_out()
            slow_tip_withdrawal(p50m, target_well, -2)
            p50m.drop_tip()
            
    
