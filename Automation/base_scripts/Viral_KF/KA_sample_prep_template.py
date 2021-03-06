import math
from opentrons.types import Point
from opentrons import protocol_api
import time
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'S2 Station A Kingfisher Version 2',
    'author': 'Hiart Maortua, Aitor Gastaminza, Arkaitz Monteju & José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Clínic Barcelona, Hospital Universitario Cruces Bilbao',
    'apiLevel': '2.0',
    'description': 'Protocol for Kingfisher sample setup viral (A)'
}

'''
'technician': $technician
'date': $date
'''

#Defined variables
##################
NUM_SAMPLES = $num_samples #excluding PC and NC
#NUM_SAMPLES = NUM_SAMPLES - 2
air_gap_vol = 0
run_id = $run_id
volume_sample = 50
lysis_volume = 100
ic_volume = 10
x_offset = [0,0]

source_type='screwcap_2ml' #'eppendorf_1.5ml' # or 'screwcap_2ml'

# Screwcap variables
diameter_screwcap = 8.25  # Diameter of the screwcap
volume_cone = 50  # Volume in ul that fit in the screwcap cone

# Calculated variables
area_section_screwcap = (math.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
screwcap_cross_section_area = math.pi * \
    diameter_screwcap**2 / 4  # screwcap cross section area


def run(ctx: protocol_api.ProtocolContext):
    import os
    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': False, 'description': 'Add lysis buffer'},
        2: {'Execute': True, 'description': 'Add samples (50ul)'},
        3: {'Execute': False, 'description': 'Add internal control (10ul)'}
    }
    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    if not ctx.is_simulating():
        # Folder and file_path for log time
        folder_path = '/var/lib/jupyter/notebooks/'+str(run_id)
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/KA_SampleSetup_viral_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none', rinse_loops = 3):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
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

    Samples = Reagent(name = 'Samples',
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 2,
                      rinse = False,
                      delay = 0,
                      reagent_reservoir_volume = 100 * 24,
                      num_wells = 24,  # num_cols comes from available columns
                      h_cono = h_cone,
                      v_fondo = volume_cone
                      )  # cone

    LBuffer = Reagent(name = 'Lysis buffer',
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      rinse = False,
                      delay = 0,
                      reagent_reservoir_volume = $Lysis_total_volume, #NUM_SAMPLES*100*1.1,
                      num_wells = $Lysis_wells,  # num_cols comes from available columns
                      h_cono = h_cone,
                      v_fondo = volume_cone
                      )  # cone

    IC = Reagent(name = 'Internal control',
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 3,
                      rinse = False,
                      delay = 0,
                      reagent_reservoir_volume = $IC_total_volume,
                      num_wells = $IC_wells,  # num_cols comes from available columns
                      h_cono = h_cone,
                      v_fondo = volume_cone
                      )  # cone

    Samples.vol_well = Samples.vol_well_original
    LBuffer.vol_well = LBuffer.vol_well_original
    IC.vol_well = IC.vol_well_original

    ##################
    # Custom functions
    def generate_source_table(source):
        '''
        Concatenate the wells from the different origin racks
        '''
        for rack_number in range(len(source)):
            if rack_number == 0:
                s = source[rack_number].wells()
            else:
                s = s + source[rack_number].wells()
        return s

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
            pipet.aspirate(post_airgap_vol, dest.top(z = 5))



    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 3, post_airgap=False, post_airgap_vol=10,
    post_dispense=False, post_dispense_vol=20,):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        if mix_height == 0:
            mix_height = 3
        pipet.aspirate(1, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        for _ in range(rounds):
            pipet.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
            pipet.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        pipet.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        if blow_out == True:
            pipet.blow_out(location.top(z=-2))  # Blow out
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, location.top(z = -2))
        if post_airgap == True:
            pipet.dispense(post_airgap_vol, location.top(z = 5))

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

    ####################################
    # load labware and modules

    ####################################
    # Load Sample racks
    if NUM_SAMPLES < 96:
        rack_num = math.ceil(NUM_SAMPLES / 24)
        ctx.comment('Used source racks are ' + str(rack_num))
        samples_last_rack = NUM_SAMPLES - rack_num * 24
    else:
        rack_num = 4

    source_tube_types={'screwcap_2ml': ['opentrons_24_tuberack_generic_2ml_screwcap','source tuberack with screwcap'],
                        'eppendorf_1.5ml': ['opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap','source tuberack with eppendorf'],
                        }

    source_racks = [ctx.load_labware(
        source_tube_types[source_type][0], slot,
        source_tube_types[source_type][1] + str(i + 1)) for i, slot in enumerate(['4', '1', '6', '3'][:rack_num])
    ]

    lysis_source_rack = ctx.load_labware('opentrons_24_tuberack_generic_2ml_screwcap', '7',
    'Lysis source')
    lysis_source = lysis_source_rack.wells()[0] # lysis comes from 1 bottle

    ic_source_rack = ctx.load_labware('opentrons_24_aluminumblock_generic_2ml_screwcap', '9',
    'Internal control source')

    ic_source = ic_source_rack.wells()[0] # internal control comes from 1 bottle

    ##################################
    # Destination plate
    dest_plate = ctx.load_labware(
        'kf_96_wellplate_2400ul', '5', 'KF 96well destination plate')

    ####################################
    # Load tip_racks
    # tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, '20µl filter tiprack')
    # for slot in ['2', '8']]
    tips300 = [ctx.load_labware('opentrons_96_filtertiprack_200ul', slot, '200µl filter tiprack')
                for slot in ['8', '11']]
    tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, '20µl filter tiprack')
                for slot in ['10']]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate

    # setup samples and destinations
    sample_sources_full = generate_source_table(source_racks)
    sample_sources = sample_sources_full[:NUM_SAMPLES]
    destinations = dest_plate.wells()[:NUM_SAMPLES]

    p20 = ctx.load_instrument(
        'p20_single_gen2', mount='left', tip_racks=tips20)
    p300 = ctx.load_instrument(
        'p300_single_gen2', mount='right', tip_racks=tips300)  # load P1000 pipette

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {p300: 0, p20: 0},  # p1000: 0},
        'maxes': {p300: len(tips300) * 96, p20: len(tips20)*96}  # ,p20: len(tips20)*96,
    }

    ############################################################################
    # STEP 1: Add lysis buffer
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        # Transfer parameters
        start = datetime.now()
        for d in destinations:
            if not p300.hw_pipette['has_tip']:
                pick_up(p300)
            # Mix the sample BEFORE dispensing
            #custom_mix(p1000, reagent = Samples, location = s, vol = volume_sample, rounds = 2, blow_out = True, mix_height = 15)
            move_vol_multichannel(p300, reagent = LBuffer, source = lysis_source, dest = d,
            vol = lysis_volume, air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 1, rinse = LBuffer.rinse, disp_height = -10,
                               blow_out = False, touch_tip = False)
            # Mix the sample AFTER dispensing
            #custom_mix(p1000, reagent = Samples, location = d, vol = volume_sample, rounds = 2, blow_out = True, mix_height = 15)
            # Drop tip and update counter
            p300.drop_tip(home_after=False)
            tip_track['counts'][p300] += 1

        # Time statistics
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] +
                    ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 2: Add Samples
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        # Transfer parameters
        start = datetime.now()
        for s, d in zip(sample_sources, destinations):
            if not p300.hw_pipette['has_tip']:
                pick_up(p300)
            # Mix the sample BEFORE dispensing
            #custom_mix(p1000, reagent = Samples, location = s, vol = volume_sample, rounds = 2, blow_out = True, mix_height = 15)
            move_vol_multichannel(p300, reagent = Samples, source = s, dest = d,
            vol=volume_sample, air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 0.4, rinse = Samples.rinse, disp_height = -40.7,
                               blow_out = False, touch_tip = True)
            # Mix the sample AFTER dispensing
            #custom_mix(p300, reagent = Samples, location = d, vol = volume_sample, rounds = 2, blow_out = True, mix_height = 15)
            # Drop tip and update counter
            p300.drop_tip(home_after=False)
            tip_track['counts'][p300] += 1

        # Time statistics
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] +
                    ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: Add internal control
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        #BEWARE, everything with the same tip!

        # Transfer parameters
        start = datetime.now()
        for d in destinations:
            if not p20.hw_pipette['has_tip']:
                pick_up(p20)
            # Mix the sample BEFORE dispensing
            #custom_mix(p1000, reagent = Samples, location = s, vol = volume_sample, rounds = 2, blow_out = True, mix_height = 15)
            move_vol_multichannel(p20, reagent = IC, source = ic_source, dest = d,
                                  vol = ic_volume, air_gap_vol = air_gap_vol, x_offset = x_offset,
                                  pickup_height = 0.4, rinse = IC.rinse, disp_height = -40.7,
                                  blow_out = True, touch_tip = False)
            # Mix the sample AFTER dispensing
            #custom_mix(p20, reagent = Samples, location = d, vol = 10, rounds = 2,
            #blow_out = True, mix_height = 2, x_offset = x_offset)
            # Drop tip and update counter
            p20.drop_tip(home_after=False)
        tip_track['counts'][p20] += 1

        # Time statistics
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] +
                    ' took ' + str(time_taken))
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

    #if not ctx.is_simulating():
        #os.system('mpg123 -f -8000 /etc/audio/speaker-test.mp3 &')
    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=False)
    ctx.home()

    ctx.comment(
        'Finished! \nMove deepwell plate (slot 2) to Station B.')

    ctx.comment('Used 200µl tips in total: ' + str(tip_track['counts'][p300]))
    ctx.comment('Used 200ul racks in total: '+str(tip_track['counts'][p300] / 96))

    ctx.comment('Used p20 tips in total: ' + str(tip_track['counts'][p20]))
    ctx.comment('Used p20 racks in total: ' + str(tip_track['counts'][p20] / 96))
