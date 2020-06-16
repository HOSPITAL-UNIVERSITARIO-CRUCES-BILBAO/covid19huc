import math
from opentrons.types import Point
from opentrons import protocol_api
import time

import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Kingfisher Pathogen Station B2 v2',
    'author': 'Aitor Gastaminza, Eva Gonzalez, José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Clínic Barcelona, Hospital Cruces Barakaldo',
    'apiLevel': '2.2',
    'description': 'Protocol for RNA extraction preparation for ThermoFisher Viral kit (ref 48383) \
    setup - sample + beads'
}

'''
'technician': '$technician',
'date': '$date'
'''

# Defined variables
##################
run_id = '43001'
NUM_SAMPLES = 48
air_gap_vol = 5

x_offset = [0,0]
L_deepwell = 8  # Deepwell side length (KingFisher deepwell)
# Screwcap variables
diameter_screwcap = 8.25  # Diameter of the screwcap
volume_cone = 50  # Volume in ul that fit in the screwcap cone

# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
screwcap_cross_section_area = math.pi * \
diameter_screwcap**2 / 4  # screwcap cross secion area
multi_well_rack_area = 8.2 * 71.2  # Cross section of the 12 well reservoir
deepwell_cross_section_area = L_deepwell**2  # deepwell cross secion area
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on


# 'kf_96_wellplate_2400ul'
def run(ctx: protocol_api.ProtocolContext):
    import os
    from opentrons.drivers.rpi_drivers import gpio
    ctx.comment('Actual used columns: ' + str(num_cols))

    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': True, 'description': 'Transfer IC'}
    }
    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    folder_path = '/var/lib/jupyter/notebooks/'+run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/Station_KB_IC_pathogen_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none',flow_rate_mix = 1, rinse_loops = 2):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_mix = flow_rate_mix
            self.rinse = bool(rinse)
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.delay = delay
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.unused=[]
            self.tip_recycling = tip_recycling
            self.vol_well_original = reagent_reservoir_volume / num_wells
            self.rinse_loops = rinse_loops


    # Reagents and their characteristics

    IC = Reagent(name='Magnetic beads and Lysis',
                    flow_rate_aspirate=1,
                    flow_rate_dispense=3,
                    rinse=True,
                    num_wells=1,
                    delay=2,
                    reagent_reservoir_volume=1050,#20 * NUM_SAMPLES * 1.1,
                    h_cono=1.95,
                    v_fondo=695)  # Prismatic


    IC.vol_well = IC.vol_well_original

    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip,
                       post_dispense=False, post_dispense_vol=20,
                       post_airgap=False, post_airgap_vol=10):
        '''
        x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        pickup_height: height from bottom where volume
        rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        blow_out, touch_tip: if True they will be done after dispensing
        '''
        # Rinse before aspirating
        if rinse == True:
            custom_mix(pipet, reagent, location = source, vol = vol,
                       rounds = 2, blow_out = True, mix_height = 0,
                       x_offset = x_offset)
        # SOURCE
        s = source.bottom(pickup_height).move(Point(x = x_offset[0]))
        pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate)  # aspirate liquid
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z = disp_height).move(Point(x = x_offset[1]))
        pipet.dispense(vol + air_gap_vol, drop,
                       rate = reagent.flow_rate_dispense)  # dispense all
        ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent
        if blow_out == True:
            pipet.blow_out(dest.top(z = -2))
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, dest.top(z = -2))
        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -5, radius = 0.9)
        if post_airgap == True:
            pipet.aspirate(post_airgap_vol, dest.top(z = -2))

    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 0.3, post_airgap=False, post_airgap_vol=10,
    post_dispense=False, post_dispense_vol=20):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        mix_trap_volumne=1
        if mix_height <= 0:
            mix_height = 0.5
        pipet.aspirate(mix_trap_volumne, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_mix)
        for _ in range(rounds):
            pipet.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_mix)
            pipet.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_mix)
        pipet.dispense(mix_trap_volumne, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_mix)
        if blow_out == True:
            pipet.blow_out(location.top(z = -2))
        if post_dispense == True:
            pipet.dispense(post_dispense_vol)
        if post_airgap == True:
            pipet.dispense(post_airgap_vol)

    def calc_height(reagent, cross_section_area, aspirate_volume, min_height = 0.3, extra_volume = 50):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume + extra_volume:
            reagent.unused.append(reagent.vol_well)
            ctx.comment('Next column should be picked')
            ctx.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
                    #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Remaining volume:' + str(reagent.vol_well))
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def distribute_custom(pipette, volume, src, dest, waste_pool, pickup_height, extra_dispensal, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        pipette.aspirate((len(dest) * volume) +
                         extra_dispensal, src.bottom(pickup_height))
        pipette.touch_tip(speed=20, v_offset=-5)
        pipette.move_to(src.top(z=5))
        pipette.aspirate(20)  # air gap
        for d in dest:
            pipette.dispense(20, d.top())
            drop = d.top(z = disp_height)
            pipette.dispense(volume, drop)
            pipette.move_to(d.top(z=5))
            pipette.aspirate(20)  # air gap

            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))

            pipette.blow_out(waste_pool.bottom(pickup_height + 3))
        return (len(dest) * volume)

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        if not ctx.is_simulating():
            if tip_track['counts'][pip] == tip_track['maxes'][pip]:
                ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0
        pip.pick_up_tip()

    def divide_destinations(l, n):
        # Divide the list of destinations in size n lists.
        for i in range(0, len(l), n):
            yield l[i:i + n]

    ####################################
    # load labware and modules
    # 12 well rack
    reagent_res = ctx.load_labware(
        'nest_12_reservoir_15ml', '5', 'Reagent deepwell plate')

    ##################################
    # Elution plate - final plate, goes to Kingfisher
    sample_plate = ctx.load_labware(
        'kf_96_wellplate_2400ul', '6',
        'KF 96 Well 2400ul elution plate')

    tips20 = [
        ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['4']
    ]

    # pipettes. P1000 currently deactivated
    m20 = ctx.load_instrument(
        'p20_multi_gen2', 'left', tip_racks=tips20)  # Load p300 multi pipette

    tip_track = {
        'counts': {m20: 0},
        'maxes': {m20: len(tips20) * 96}
    }

    # Divide destination wells in small groups for P300 pipette
    #destinations = list(divide_destinations(sample_plate.wells()[:NUM_SAMPLES], size_transfer))
    IC.reagent_reservoir = reagent_res.rows()[0][10]  # 1 row, 4 columns (first ones)
    #pipette multiple  definition jump example play
    work_destinations_cols = sample_plate.rows()[0][:num_cols]


    ############################################################################
    # STEP 1: TRANSFER IC
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        # Transfer parameters
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        IC_transfer_vol = [10]  # Two rounds of 130
        rinse = True
        for i in range(num_cols):
            if not m20.hw_pipette['has_tip']:
                pick_up(m20)
            for j, transfer_vol in enumerate(IC_transfer_vol):
                # Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(
                    reagent = IC, cross_section_area = multi_well_rack_area,
                    aspirate_volume = transfer_vol * 8, min_height=0.2, extra_volume=5)

                ctx.comment(
                    'Aspirate from reservoir column: ' + str(IC.col))
                ctx.comment('Pickup height is ' + str(0.2))

                rinse = False

                move_vol_multichannel(m20, reagent=IC, source=IC.reagent_reservoir,
                                      dest=work_destinations_cols[i], vol=transfer_vol,
                                      air_gap_vol=air_gap_vol, x_offset=x_offset,
                                      pickup_height=0.2, disp_height = -40.7,
                                      rinse=rinse, blow_out = True, touch_tip=False, post_airgap=True)


                m20.drop_tip(home_after=False)
                tip_track['counts'][m20] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # Export the time log to a tsv file
    if not ctx.is_simulating():
        with open(file_path, 'w') as f:
            f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
            for key in STEPS.keys():
                row = str(key)
                for key2 in STEPS[key].keys():
                    row += '\t' + format(STEPS[key][key2])
                f.write(row + '\n')
        f.close()


    ############################################################################
    # Light flash end of program

    import os
    if not ctx.is_simulating():
        os.system('mpg123 -f -8000 /etc/audio/speaker-test.mp3 &')
    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        #ctx._hw_manager.hardware.set_button_light(1,0,0)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        #ctx._hw_manager.hardware.set_button_light(0,0,1)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=False)
    #ctx._hw_manager.hardware.set_button_light(0,1,0)

    ctx.comment('Finished! \nMove plate to KingFisher')
