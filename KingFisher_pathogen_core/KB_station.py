import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import os
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Kingfisher Viral Station',
    'author': 'Aitor Gastaminza,  José Luis Villanueva & Eva González (jlvillanueva@clinic.cat)',
'source': 'Hospital Clínic Barcelona, Hospital Cruces Bilbao',
    'apiLevel': '2.0',
    'description': 'Protocol to fill KingFisher Deepwell plates with reagents - Pathogen Kit (ref 4462359)'
}

'''
'technician': '$technician'
'date': '$date'
'''


#Defined variables
##################

NUM_SAMPLES = 96
air_gap_vol = 10
air_gap_vol_elutionbuffer = 10
ic_volume = 10
lysis_vol = 260
beads_vol = 260
elution_buffer_vol = 90
wash_buffer1_vol = 300
wash_buffer2_vol = 300
waiting = 10 # minutes
pipette_allowed_capacity=180
max_multiwell_volume = 13500


run_id =  '$run_id'

x_offset = [0,0]
multi_well_rack_area = 8.2 * 71.2  # Cross section of the 12 well reservoir
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on


def run(ctx: protocol_api.ProtocolContext):
    ctx.comment('Actual used columns: ' + str(num_cols))
    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': True, 'description': 'Add 10 ul IC'},
        2: {'Execute': True, 'description': 'Add 260 ul Lysis Buffer'},
        3: {'Execute': True, 'description': 'Wait for 10 minutes', 'wait_time': 600}, # 10 minutes of waiting
        4: {'Execute': True, 'description': 'Add 260 ul Beads'},
        5: {'Execute': True, 'description': 'Add 300 ul Wash Buffer 1 - Round 1'},
        6: {'Execute': True, 'description': 'Add 450 ul Wash Buffer 2 - Round 1'},
        7: {'Execute': True, 'description': 'Add 90 ul Elution Buffer'}
        }

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0
    folder_path = '/var/lib/jupyter/notebooks/'+run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/KB_PlateFilling_pathogen_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none', rinse_loops = 2):

            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.rinse = bool(rinse)
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.delay = delay #Delay of reagent in dispense
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
    WashBuffer1 = Reagent(name='Wash Buffer 1',
                          flow_rate_aspirate=0.75,
                          flow_rate_dispense=1,
                          rinse=True,
                          delay=2,
                          reagent_reservoir_volume=wash_buffer2_vol*NUM_SAMPLES*1.1,
                          num_wells=math.ceil((NUM_SAMPLES + 5) * wash_buffer2_vol / max_multiwell_volume),
                          h_cono=1.95,
                          v_fondo=695)  # Flat surface

    WashBuffer2 = Reagent(name='Wash Buffer 2',
                          flow_rate_aspirate=0.75,
                          flow_rate_dispense=1,
                          rinse=True,
                          delay=2,
                          reagent_reservoir_volume=wash_buffer1_vol*1.1*NUM_SAMPLES,
                          num_wells=math.ceil((NUM_SAMPLES + 5) * wash_buffer1_vol / max_multiwell_volume),
                          h_cono=1.95,
                          v_fondo=695)  # Flat surface

    Lysis = Reagent(name='Lysis Buffer',
                          flow_rate_aspirate=0.75,
                          flow_rate_dispense=0.5,
                          rinse=False,
                          delay=2,
                          reagent_reservoir_volume=lysis_vol*1.1*NUM_SAMPLES,
                          num_wells=math.ceil((NUM_SAMPLES + 5) * lysis_vol / max_multiwell_volume),
                          h_cono=1.95,
                          v_fondo=695)  # Flat surface

    ElutionBuffer = Reagent(name='Elution Buffer',
                            flow_rate_aspirate=1,
                            flow_rate_dispense=1,
                            rinse=False,
                            delay=0,
                            reagent_reservoir_volume=elution_buffer_vol*1.1*NUM_SAMPLES,#50*NUM_SAMPLES,
                            num_wells=math.ceil((NUM_SAMPLES + 5) * elution_buffer_vol / max_multiwell_volume),
                            h_cono=1.95,
                            v_fondo=695)  # Prismatic

    IC = Reagent(name='Internal control',
                    flow_rate_aspirate=1,
                    flow_rate_dispense=3,
                    rinse=True,
                    num_wells=1,
                    delay=2,
                    reagent_reservoir_volume=1800,#20 * NUM_SAMPLES * 1.1,
                    h_cono=1.95,
                    v_fondo=695)  # Prismatic

    Beads = Reagent(name='Magnetic beads',
                    flow_rate_aspirate=0.5,
                    flow_rate_dispense=0.5,
                    rinse=True,
                    num_wells=math.ceil((NUM_SAMPLES + 5) * beads_vol / max_multiwell_volume),
                    delay=2,
                    reagent_reservoir_volume=beads_vol*1.15*NUM_SAMPLES,#20 * NUM_SAMPLES * 1.1,
                    h_cono=1.95,
                    v_fondo=695)  # Prismatic


    Beads.vol_well = Beads.vol_well_original
    IC.vol_well = IC.vol_well_original
    WashBuffer1.vol_well = WashBuffer1.vol_well_original
    WashBuffer2.vol_well = WashBuffer2.vol_well_original
    ElutionBuffer.vol_well = ElutionBuffer.vol_well_original
    Lysis.vol_well = Lysis.vol_well_original
    ctx.comment(str(Beads.vol_well))

    ##################
    # Custom functions
    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip,
                       post_dispense=False, post_dispense_vol=20,
                       post_airgap=True, post_airgap_vol=10):
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
                       rounds = reagent.rinse_loops, blow_out = True, mix_height = 0,
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
        if post_airgap == True:
            pipet.dispense(post_airgap_vol, dest.top(z = -2))
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, dest.top(z = -2))
        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -5, radius = 0.9)


    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 3, post_airgap=True, post_airgap_vol=10,
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
            ctx.comment('Remaining volume now will be:' + str(reagent.vol_well))
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

    def divide_volume(volume,max_vol):
        num_transfers=math.ceil(volume/max_vol)
        vol_roundup=math.ceil(volume/num_transfers)
        last_vol=volume-vol_roundup*(num_transfers-1)
        vol_list=[vol_roundup for v in range(1,num_transfers)]
        vol_list.append(last_vol)
        return vol_list

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
    ##########

    def find_side(col):
        '''
        Detects if the current column has the magnet at its left or right side
        '''
        if col % 2 == 0:
            side = -1  # left
        else:
            side = 1
        return side


####################################
    # load labware and modules

    # 12 well rack
    ####################################
    reagent_res = ctx.load_labware(
        'nest_12_reservoir_15ml', '1', 'Reservoir 12 channel, column 1')

    ic_res = ctx.load_labware(
        'nest_12_reservoir_15ml', '4', 'Reservoir 12 channel, column 1')


    # Wash Buffer 1 100ul Deepwell plate
    ############################################
    WashBuffer1_100ul_plate1 = ctx.load_labware(
        'kf_96_wellplate_2400ul', '9', 'Wash Buffer 1 Deepwell plate 1')

    # Wash Buffer 2 100ul Deepwell plate
    ############################################
    WashBuffer2_100ul_plate1 = ctx.load_labware(
        'kf_96_wellplate_2400ul', '6', 'Wash Buffer 2 Deepwell plate 1')

    # Plate with samples
    ############################################
    kf_plate = ctx.load_labware(
        'kf_96_wellplate_2400ul', '2', 'Deepwell plate with samples')

    # Elution Deepwell plate
    ############################################
    ElutionBuffer_50ul_plate = ctx.load_labware(
        'kingfisher_std_96_wellplate_550ul', '10', 'Elution Buffer 90 ul STD plate')


####################################
    # Load tip_racks
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
               for slot in ['7']]
    tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, '20µl filter tiprack')
               for slot in ['8']]

################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    WashBuffer1.reagent_reservoir = reagent_res.rows(
    )[0][:3] # position 1, 3 columns

    WashBuffer2.reagent_reservoir = reagent_res.rows(
    )[0][3:7] #position 2, 3 columns

    Lysis.reagent_reservoir = reagent_res.rows(
    )[0][7:9] #position 3, 2 columns

    Beads.reagent_reservoir = reagent_res.rows(
    )[0][9:11] #position 4, 2 columns

    ElutionBuffer.reagent_reservoir = reagent_res.rows(
    )[0][11:] #position 5, 1 column

    IC.reagent_reservoir = ic_res.rows(
    )[0][:1] #position 1

    # columns in destination plates to be filled depending the number of samples
    wb1plate1_destination = WashBuffer1_100ul_plate1.rows()[0][:num_cols]
    wb2plate1_destination = WashBuffer2_100ul_plate1.rows()[0][:num_cols]
    elutionbuffer_destination = ElutionBuffer_50ul_plate.rows()[0][:num_cols]
    kf_destination = kf_plate.rows()[0][:num_cols]

    # pipette
    m300 = ctx.load_instrument(
        'p300_multi_gen2', 'right', tip_racks=tips300)  # Load multi pipette

    m20 = ctx.load_instrument(
        'p20_multi_gen2', 'left', tip_racks=tips20)  # Load multi pipette

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0, m20: 0},
        'maxes': {m300: len(tips300)*96, m20: len(tips20)*96}
    }

    ############################################################################
    # STEP 1 Filling with IC
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        rinse = False  # Only first time
        ########
        # Wash buffer dispense
        for i in range(num_cols):
            if not m20.hw_pipette['has_tip']:
                pick_up(m20)
            move_vol_multichannel(m20, reagent = IC, source = IC.reagent_reservoir[IC.col],
                           dest = kf_destination[i], vol = ic_volume,
                           air_gap_vol = air_gap_vol, x_offset = x_offset,
                           pickup_height = 1, rinse = IC.rinse, disp_height = -41,
                           blow_out = False, touch_tip = False, post_airgap=True)
            m20.drop_tip(home_after=False)
            tip_track['counts'][m20] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 2 Add Lysis
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        rinse = False  # Only first time

        if (lysis_vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
        # calculate what volume should be transferred in each step
            vol_list=divide_volume(lysis_vol, pipette_allowed_capacity)

        ########
        # Wash buffer dispense
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j, transfer_vol in enumerate(vol_list):
                '''if (i == 0 and j == 0):
                    rinse = True
                else:
                    rinse = False'''
                [pickup_height,col_change]=calc_height(Lysis, multi_well_rack_area, transfer_vol*8)
                move_vol_multichannel(m300, reagent = Lysis, source = Lysis.reagent_reservoir[Lysis.col],
                               dest = kf_destination[i], vol = transfer_vol,
                               air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 1, rinse = Lysis.rinse, disp_height = -2,
                               blow_out = True, touch_tip = False, post_airgap=True)

        m300.drop_tip(home_after=False)
        tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3 Wait time
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer magnetic beads
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Waiting for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # minutes=2
        ctx.comment(' ')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)

    ############################################################################
    # STEP 4 Filling with Beads
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        if (beads_vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
        # calculate what volume should be transferred in each step
            vol_list=divide_volume(beads_vol, pipette_allowed_capacity)

        ########
        # Wash buffer dispense
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j, transfer_vol in enumerate(vol_list):
                if (i == 0 and j == 0):
                    rinse = True
                else:
                    rinse = False
                [pickup_height,col_change]=calc_height(Beads, multi_well_rack_area, transfer_vol*8)
                move_vol_multichannel(m300, reagent = Beads, source = Beads.reagent_reservoir[Beads.col],
                               dest = kf_destination[i], vol = transfer_vol,
                               air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = pickup_height, rinse = rinse, disp_height = -2,
                               blow_out = True, touch_tip = False, post_airgap=True)

        m300.drop_tip(home_after=False)
        tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 5  Filling with WashBuffer1
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        if (wash_buffer1_vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
        # calculate what volume should be transferred in each step
            vol_list=divide_volume(wash_buffer1_vol, pipette_allowed_capacity)

        ########
        # Wash buffer dispense
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j, transfer_vol in enumerate(vol_list):
                if (i == 0 and j == 0):
                    rinse = True
                else:
                    rinse = False
                [pickup_height,col_change]=calc_height(WashBuffer1, multi_well_rack_area, transfer_vol*8)
                move_vol_multichannel(m300, reagent = WashBuffer1, source = WashBuffer1.reagent_reservoir[WashBuffer1.col],
                               dest = wb1plate1_destination[i], vol = transfer_vol,
                               air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 1, rinse = rinse, disp_height = -2,
                               blow_out = True, touch_tip = False, post_airgap=True)

        m300.drop_tip(home_after=False)
        tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 6  Filling with WashBuffer2
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        if (wash_buffer2_vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
        # calculate what volume should be transferred in each step
            vol_list=divide_volume(wash_buffer2_vol, pipette_allowed_capacity)

        ########
        # Wash buffer dispense
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j, transfer_vol in enumerate(vol_list):
                if (i == 0 and j == 0):
                    rinse = True
                else:
                    rinse = False
                [pickup_height,col_change]=calc_height(WashBuffer2, multi_well_rack_area, transfer_vol*8)
                move_vol_multichannel(m300, reagent = WashBuffer2, source = WashBuffer2.reagent_reservoir[WashBuffer2.col],
                               dest = wb2plate1_destination[i], vol = transfer_vol,
                               air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 1, rinse = rinse, disp_height = -2,
                               blow_out = True, touch_tip = False, post_airgap=True)

        m300.drop_tip(home_after=False)
        tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)



    ############################################################################
    # STEP 7 Transfer Elution buffer
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        # Elution buffer

        ########
        # Water or elution buffer
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
                # Calculate pickup_height based on remaining volume and shape of container
            ctx.comment(
                'Aspirate from Reservoir column: ' + str(ElutionBuffer.col))
            move_vol_multichannel(m300, reagent = ElutionBuffer, source = ElutionBuffer.reagent_reservoir[ElutionBuffer.col],
                          dest = elutionbuffer_destination[i], vol = elution_buffer_vol,
                          air_gap_vol = air_gap_vol_elutionbuffer, x_offset = x_offset,
                          pickup_height = 0.5, rinse = False, disp_height = -4,
                          blow_out = True, touch_tip = False, post_airgap=True)
        m300.drop_tip(home_after=False)
        tip_track['counts'][m300] += 8
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
    if not ctx.is_simulating():
        from opentrons.drivers.rpi_drivers import gpio

        for i in range(3):
            ctx._hw_manager.hardware.set_lights(rails=False)
            ctx._hw_manager.hardware.set_lights(button=(1,0,0))
            time.sleep(0.3)
            ctx._hw_manager.hardware.set_lights(rails=True)
            ctx._hw_manager.hardware.set_lights(button=(0,0,1))
            time.sleep(0.3)
            ctx._hw_manager.hardware.set_lights(rails=False)
        ctx._hw_manager.hardware.set_lights(button=(0,1,0))

        ctx.comment(
            'Finished! \nMove deepwell plates to KingFisher extractor.')
        ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
        ctx.comment('Used racks in total: ' + str(tip_track['counts'][m300] / 96))
