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
    NUM_SAMPLES = 8  # Define the number of samples (up to 48)
    NUM_COLUMNS = (NUM_SAMPLES + 7) // 8  # Calculate number of columns for 96-well plates
    
    # TESTING PARAMETERS
    DRY_RUN = 0.01      # use 0.01 to shorten wait times if it is dry run, otherwise 1
    TIPRECYCLE = True   # change to False if not a dry run (eg don't recycle tips)

    skip_mixbeadsandrna = False     # Toggle when testing certain blocks, same below
    skip_thermocycler1 = True
    skip_beadmix = True
    skip_supremoval = False
    skip_washes = False
    skip_1stelution = True
    skip_thermocycler2 = True
    skip_rebinding = True
    skip_supremoval2 = True
    skip_finalwash = True
    skip_washremoval = True
    skip_2ndelution = True


    # EXPERIMENTAL PARAMETERS
    deadvol_reservoir = 2000
    deadvol_plate = 5
    clearance_reservoir = 2
    clearance_bead_pellet = 2
    clearance_beadresuspension = 3
    offset_x = 1

    #======== DECK SETUP ========
    # COLUMN 1 
    thermocycler = protocol.load_module('thermocycler module gen2')
    temp_block = protocol.load_module('temperature module gen2', 'D1')
    temp_adapter = temp_block.load_adapter('opentrons_96_well_aluminum_block')     
    reagent_plate = temp_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    
    # COLUMN 2 
    mag_block = protocol.load_module('magneticBlockV1', 'A2')
    tiprack_50ul = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B2')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    sample_plate = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", "D2")
     
    # COLUMN 3 
    tiprack_200ul_1 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'A3')
    tiprack_200ul_2 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'B3')
    if TIPRECYCLE:
        waste = protocol.load_trash_bin("C3")
    else:
        waste_chute = protocol.load_waste_chute()
     
    # COLUMN 4 
    # empty
 
    #======== PIPETTES ========
    p1000m = protocol.load_instrument('flex_8channel_1000', 'left', tip_racks=[tiprack_200ul_1,tiprack_200ul_2])
    p50m = protocol.load_instrument('flex_8channel_50', 'right', tip_racks=[tiprack_50ul])
     
    #======== REAGENTS ========
    # reagent plate
    # binding buffer
    rna_binding_buffer = [
     reagent_plate.columns_by_name()[name] for name in ['1', '2']]  # RNA Binding Buffer (on ice)
    
    # binding buffer (2X 50 uL per sample) - i
    for index, column in enumerate(rna_binding_buffer):
        if not index:
            num = NUM_COLUMNS if NUM_COLUMNS <= 3 else 3    # only 1 column for up to 3 samples
        else:
            num = NUM_COLUMNS - 3 if NUM_COLUMNS > 3 else 0 # 2 columns for more than 3 samples
        # sets volume needed in each column
        if num:
            column[0].liq_vol = 50*num + deadvol_plate
        else:
            column[0].liq_vol = 0

    # reagent reservoir
    ### TODO check if volume correct - wash buffer 180 uL per sample - for three washes
    wash_buffer = [reservoir.wells_by_name()[well] for well in ['A1','A2']]
    for well in wash_buffer:
        well.liq_vol = 180*NUM_COLUMNS*p1000m.channels + deadvol_reservoir
    
    [tris] = [reservoir.wells_by_name()[well] for well in ['A4']]
    waste = [reservoir.wells_by_name()[well] for well in ['A9', 'A10', 'A11', 'A12']]

    #======== DEFINING LIQUIDS ========   
    sample_wells = protocol.define_liquid(name="Samples",
                                          description="50 µL sample + 50 µL 2X binding buffer",
                                          display_color="#C0C0C0") # Silver
    sample_vol = 100    # 50 µL RNA sample with 50 µL 2x binding buffer added just prior

    wash_buffer_wells = protocol.define_liquid(name="wash buffer", 
                                               description="wash buffer, 180µl * sample * number of washes", 
                                               display_color="#0000FF") # Blue
    wash_buffer_vol = 180 * NUM_SAMPLES * 3

    tris_wells = protocol.define_liquid(name="Tris buffer",
                                        description="50 µL per sample * sample count",
                                        display_color="#FFA500") # Orange
    tris_vol = 50 * NUM_SAMPLES

    binding_buffer_wells = protocol.define_liquid(name="2X RNA binding buffer",
                                        description="50 µL per sample * sample count",
                                        display_color="#008000") # Green
    binding_buffer_vol = 50 * NUM_SAMPLES

    #======== LOADING LIQUIDS ========
    for well in sample_plate.wells()[:NUM_SAMPLES]:
        well.load_liquid(liquid=sample_wells, volume=sample_vol)

    for column in reservoir.columns()[0:2]:
        for well in column: 
            well.load_liquid(liquid=wash_buffer_wells, volume=wash_buffer_vol)

    for well in reservoir.columns()[3]:
        well.load_liquid(liquid=tris_wells, volume=tris_vol)

    for well in reagent_plate.columns()[0]:
        well.load_liquid(liquid=binding_buffer_wells, volume=binding_buffer_vol)

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

    # FUNCTION 3: BEAD MIXING
    def bead_mixing(reps, vol, aspirate_rate=0.8, dispense_rate=2):
        for index, column in enumerate(sample_plate.columns()[:NUM_COLUMNS]):
            pick_up_or_refill(p1000m)
            
            # number of mixing reps
            for rep in range(reps):
                p1000m.aspirate(vol, column[0].bottom(clearance_bead_pellet), rate=aspirate_rate)
                p1000m.dispense(vol, column[0].bottom(clearance_beadresuspension), rate=dispense_rate)

                if rep == reps - 1:
                    # side touch with blowout after last mix
                    p1000m.move_to(
                        column[0].top(-2).move(types.Point(
                        x=column[0].diameter / 2, y=0, z=0))) # TODO try speed=0.5
                    p1000m.blow_out()
                    protocol.delay(seconds=1.5)
                    p1000m.move_to(column[0].top())
            p1000m.drop_tip()

    # FUNCTION 4: REMOVE SUPERNATANT
    # TODO Double check the top(), bottom(2) and rate parameters
    def remove_sup(volume1,volume2,waste_loc):       
        pick_up_or_refill(p1000m)
        
        p1000m.move_to(column[0].top())
        p1000m.air_gap(20)
        p1000m.aspirate(volume1, column[0].bottom(4), rate=0.2)
        p1000m.aspirate(volume2, 
                        column[0].bottom(clearance_bead_pellet).move(types.Point(x=offset_x, y=0, z=0)),
                        rate=0.05)
        protocol.delay(seconds=1)

        p1000m.dispense(volume1+volume2+20, waste[waste_loc].top(), rate=2)
        protocol.delay(seconds=1)
        p1000m.blow_out()
        p1000m.air_gap(20)
        p1000m.drop_tip()

    # FUNCTION 5: WASH MIXING
    def wash_mixing(volume,reps):
        ### old mixing protocol with alternate dispense location
        # for rep in range(reps):
                    
            # clearance_mixdispense = 6 if (rep % 2) else clearance_beadresuspension
            # offset_x_mixdispense = -1 if rep % 2 else offset_x
            # loc = column[0].bottom(clearance_mixdispense).move(types.Point(x=offset_x_mixdispense, y=0, z=0))
    
            # # aspirate and dispense at different locations
            # p1000m.aspirate(volume, column[0].bottom(1))
            # p1000m.dispense(volume, loc, rate=3)

        p1000m.mix(reps, volume, column[0].bottom(3), rate=0.4)
        slow_tip_withdrawal(p1000m, column[0])
        p1000m.blow_out(column[0].bottom(13))
        p1000m.touch_tip(radius=0.8, v_offset=-3, speed=5)
        p1000m.drop_tip()

    
    # FUNCTION X: LIQUID HEIGHT IN A WELL
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
    protocol.delay(minutes=2*DRY_RUN)

    # STEP 19: Remove and discard supernatant
    protocol.comment("Step 19: Remove and discard supernatant.")

    if not skip_supremoval:
        for index, column in enumerate(sample_plate.columns()[:NUM_COLUMNS]):
            remove_sup(75,75,waste_loc=-2)
    
    # STEP 20-25: Add wash buffer and wash twice
    protocol.comment("Step 20-25: Take off magnet, wash beads and mix 10 times. Repeat")

    if not skip_washes:

        # washing twice total
        for repeat in range(2):
            
            # move plate off magnet for washes
            protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

            # dispense wash buffer to all bead pellets
            pick_up_or_refill(p1000m)
            for column in sample_plate.columns()[:NUM_COLUMNS]:
                p1000m.aspirate(180, wash_buffer[repeat].bottom(clearance_reservoir)) # ht instead of 2 if using
                p1000m.dispense(180, column[0].top(), rate = 0.8)
                # protocol.delay(seconds=1)
                # p1000m.blow_out()

            # mix all
            for column in sample_plate.columns()[:NUM_COLUMNS]:
                if not p1000m.has_tip:
                    pick_up_or_refill(p1000m)
                wash_mixing(volume=144,reps=15)

            # move plate to magnet for supernatant removal
            protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
            protocol.delay(minutes=2*DRY_RUN)

            # remove supernatant
            for column in sample_plate.columns()[:NUM_COLUMNS]:
                remove_sup(130,50,repeat)
                
            # complete removal of last wash
            if repeat:
                for column in sample_plate.columns()[:NUM_COLUMNS]:
                    
                    pick_up_or_refill(p1000m)
                    
                    # loop to move closer and closer to the bottom
                    # TODO check that the clearance of 0.7-0 isn't too low??
                    for clearance in [0.7, 0.4, 0.2, 0]:
                        loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0,z=0))
                        p1000m.aspirate(25, loc)
                    p1000m.drop_tip()
            
    protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 26: adding Tris for elution
    protocol.comment("Step 26: First elution with Tris, mix 10 times.")

    if not skip_1stelution:

        for column in sample_plate.columns()[:NUM_COLUMNS]:
            
            pick_up_or_refill(p1000m)
            p1000m.aspirate(50, tris.bottom(clearance_reservoir))
            
            # dispense Tris quickly on the side of pellet
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))            
            p1000m.dispense(50, loc, rate=2.5)
            
            # Mixing loop - alternating vertical clearance and horizontal offset
            wash_mixing(volume=43,reps=15)

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
    
    # STEP 29-32: Add 50 μl of RNA Binding Buffer (2X), mix and re-bind, repeat 1x
    
    protocol.comment("Step 29-32: Add RNA binding buffer and mix, re-bind, repeat 1x")
    
    if not skip_rebinding:
        source = rna_binding_buffer[0][0]

        for index, column in enumerate(sample_plate.columns()[:NUM_COLUMNS]):
            pick_up_or_refill(p1000m)

            # use second column of Binding buffer for sample columns >3
            if index > 2:
                source = rna_binding_buffer[1][0]
            
            p1000m.aspirate(50, source.bottom(1))
            p1000m.dispense(50, column[0].bottom(2))

            # mixing
            for rep in range(10):
                p1000m.aspirate(80, column[0].bottom(2))
                p1000m.dispense(80, column[0].bottom(2))
                
                # last mixing rep - slow tip withdrawal, move slightly offcenter at top of well (top(-2)
                # blow_out() to ensure no remaining droplets in tip
                if rep == 9:
                    protocol.delay(seconds=1)
                    p1000m.move_to(well.top(),speed = 5)

                    p1000m.move_to(
                     column[0].top(-2).move(types.Point(
                      x=column[0].diameter / 2, y=0, z=0)))
                    p1000m.blow_out()
                    p1000m.move_to(column[0].top())
            p1000m.drop_tip()

        protocol.delay(minutes=5*DRY_RUN)
        bead_mixing(10, 80, aspirate_rate=0.2, dispense_rate=0.2) # 10 mixes with 80µl
        protocol.delay(minutes=5*DRY_RUN)
    
    # STEP 33: move plate to magnet for 2 min to separate the poly(A) RNA
    protocol.comment("Step 33: Second binding - separate RNA from solution.")
    
    protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=2*DRY_RUN)

    # STEP 34-35: Remove and discard supernatant, move off magnet
    protocol.comment("Step 34: Remove and discard supernatant.")
    if not skip_supremoval2:
        for column in sample_plate.columns()[:NUM_COLUMNS]:
            remove_sup(75,75,waste_loc=-2)

    protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 36: Add final wash buffer and mix
    protocol.comment("Step 36: Wash the beads once with 200 μl of Wash Buffer, mix thoroughly.")

    if not skip_finalwash:
        for index, column in enumerate(sample_plate.columns()[:NUM_COLUMNS]):

            pick_up_or_refill(p1000m)

            ### TODO this section tracks the liquid height in the reservoir, check if needd 
            # wash_buffer[-1].liq_vol -= 180*p1000m.channels
            # ht = liq_height(
            #   wash_buffer[-1]) - 3 if liq_height(
            #   wash_buffer[-1]) - 3 > 1 else 1
            p1000m.aspirate(180, wash_buffer[-1].bottom(2))
            p1000m.dispense(180, column[0].bottom(2)) # ht instead of 2 if using
            wash_mixing(volume=144,reps=15)

    # STEP 37: Place tubes on magnetic rack
    protocol.comment("STEP 37: Place the tubes on the magnetic rack for 2 minutes")
    protocol.move_labware(labware=sample_plate, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=2*DRY_RUN)

    # STEP 38-39: Remove wash buffer completely and move plate off magnet
    protocol.comment("STEP 38-39: Remove and discard wash buffer, remove from magnet")

    # TODO Double check the top(), bottom(2) and rate parameters
    if not skip_washremoval:
        for column in sample_plate.columns()[:NUM_COLUMNS]:
            pick_up_or_refill(p1000m)
            
            p1000m.move_to(column[0].top())
            p1000m.air_gap(20)
            p1000m.aspirate(110, column[0].bottom(4), rate=0.2)
            
            p1000m.aspirate(
                70, column[0].bottom(1).move(types.Point(x=offset_x, y=0, z=0)), rate=0.2)
            
            p1000m.dispense(200, waste[-1].top(), rate=2)
            protocol.delay(seconds=1)
            p1000m.blow_out()

            p1000m.drop_tip()
        
        # complete removal of liquid
        for column in sample_plate.columns()[:NUM_COLUMNS]:

            pick_up_or_refill(p1000m)
            
            # loop to move closer to the bottom
            # TODO check that 0 is not too close
            for clearance in [0.7, 0.4, 0.2, 0]:
                loc = column[0].bottom(clearance).move(types.Point(x=offset_x, y=0, z=0))
                p1000m.aspirate(25, loc)
            p1000m.drop_tip()

    protocol.move_labware(labware=sample_plate, new_location="D2", use_gripper=True)

    # STEP 40: Last elution
    ## TODO This might need to be skipped to be done by hand
    if not skip_2ndelution:
        for index, column in enumerate(sample_plate.columns()[:NUM_COLUMNS]):
            
            pick_up_or_refill(p50m)
            p50m.aspirate(11.5, tris.bottom(clearance_reservoir))
            
            # calculate dispense location with alternating offset
            loc = column[0].bottom(clearance_beadresuspension).move(types.Point(x=offset_x, y=0, z=0))
            
            # dispense Tris carefully
            p50m.dispense(11.5, loc)
            protocol.delay(seconds=1)
            
            p50m.move_to(well.top(),speed = 5)
            p50m.move_to(
             column[0].top(-2).move(types.Point(
              x=column[0].diameter / 2, y=0, z=0)))
            p50m.blow_out()
            p50m.move_to(column[0].top())
            p50m.drop_tip()


            