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
    'protocolName': 'B Station',
    'author': 'Hiart Maortua & Aitor Gastaminza based on codes developed by José Luis Villanueva & Eva González (jlvillanueva@clinic.cat)',
'source': 'Hospital Clínic Barcelona, Hospital Cruces Bilbao',
    'apiLevel': '2.0',
    'description': 'Protocol for RNA extraction using MagMax Viral reactives'
}

'''
'technician': '$technician'
'date': '$date'
'''


#Defined variables
##################

NUM_SAMPLES = 24
air_gap_vol = 10
air_gap_vol_elutionbuffer = 10
run_id =  '$run_id'
air_gap_ic = 5
WBone_vol=100
WBtwo_vol=100
elution_vol=50
ic_vol = 10
beads_vol = 20
mag_height = 6.1

x_offset_rs = 1 #Offset of the pickup when magnet is ON
multi_well_rack_area = 8.2 * 71.2  # Cross section of the 12 well reservoir
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on


def run(ctx: protocol_api.ProtocolContext):
    x_offset=[0,0]
    ctx.comment('Actual used columns: ' + str(num_cols))
    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': True, 'description': 'Add 100 ul Lysis Buffer and then move to station A'},
        2: {'Execute': True, 'description': 'Add IC to the plate coming from station A'},
        3: {'Execute': True, 'description': 'Mix beads'},
        4: {'Execute': True, 'description': 'Transfer beads and mix (dispose tip)'},
        5: {'Execute': True, 'description': 'Wait with magnet OFF after beads', 'wait_time': 10},  # 60
        6: {'Execute': True, 'description': 'Wait with magnet ON after beads', 'wait_time': 10},  # 900
        7: {'Execute': True, 'description': 'Remove supernatant'},
        8: {'Execute': True, 'description': 'Add W1 and mix (magnet OFF)'},
        9: {'Execute': True, 'description': 'Wait with magnet ON', 'wait_time': 900},  # 900
        10: {'Execute': True, 'description': 'Remove W1'},
        11: {'Execute': True, 'description': 'Add W2 and mix (magnet OFF)'},
        12: {'Execute': True, 'description': 'Wait with magnet ON', 'wait_time': 900},  # 900
        13: {'Execute': True, 'description': 'Remove supernatant'},
        14: {'Execute': True, 'description': 'Allow to dry (magnet ON)', 'wait_time': 300},
        15: {'Execute': True, 'description': 'Add 50 ul Elution Buffer (magnet OFF) and mix'},
        16: {'Execute': True, 'description': 'Wait with magnet OFF after water', 'wait_time': 60},  # 60
        17: {'Execute': True, 'description': 'Wait with magnet ON after water', 'wait_time': 120},  # 300
        18: {'Execute': True, 'description': 'Transfer to final elution plate'}

        }

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0
    folder_path = '/var/lib/jupyter/notebooks/'+run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/B_station_manual_viral_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none', rinse_loops = 3, flow_rate_dispense_mix = 3, flow_rate_aspirate_mix = 3):

            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_dispense_mix
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
                          reagent_reservoir_volume=3500, #100*NUM_SAMPLES,
                          num_wells=1,
                          h_cono=1.95,
                          v_fondo=695,
                          tip_recycling = ['4'])  # Flat surface

    WashBuffer2 = Reagent(name='Wash Buffer 2',
                          flow_rate_aspirate=0.75,
                          flow_rate_dispense=1,
                          rinse=True,
                          delay=2,
                          reagent_reservoir_volume=3500, #100*NUM_SAMPLES,
                          num_wells=1,
                          h_cono=1.95,
                          v_fondo=695,
                          tip_recycling = ['7'])  # Flat surface

    Lysis = Reagent(name='Lysis Buffer',
                          flow_rate_aspirate=0.75,
                          flow_rate_dispense=0.5,
                          rinse=False,
                          rinse_loops=3,
                          delay=2,
                          reagent_reservoir_volume=3500, #100*NUM_SAMPLES,
                          num_wells=1,
                          h_cono=1.95,
                          v_fondo=695,
                          tip_recycling = ['2'])  # Flat surface

    ElutionBuffer = Reagent(name='Elution Buffer',
                            flow_rate_aspirate=1,
                            flow_rate_dispense=1,
                            rinse=False,
                            delay=0,
                            reagent_reservoir_volume=1600,#50*NUM_SAMPLES,
                            num_wells=1,
                            h_cono=1.95,
                            v_fondo=695)  # Prismatic

    Elution = Reagent(name='Elution',
                      flow_rate_aspirate=1,
                      flow_rate_dispense=1,
                      rinse=False,
                      reagent_reservoir_volume=50,
                      delay=0,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=4,
                      v_fondo=4 * math.pi * 8**2/4 /3)  # Sphere

    IC = Reagent(name='Magnetic beads and Lysis',
                    flow_rate_aspirate=1,
                    flow_rate_dispense=3,
                    rinse=False,
                    num_wells=1,
                    delay=2,
                    reagent_reservoir_volume=1050,#20 * NUM_SAMPLES * 1.1,
                    h_cono=1.95,
                    v_fondo=695)  # Prismatic

    Beads = Reagent(name='Magnetic beads and Lysis',
                    flow_rate_aspirate=0.5,
                    flow_rate_dispense=0.5,
                    flow_rate_aspirate_mix=4,
                    flow_rate_dispense_mix=6,
                    rinse=True,
                    rinse_loops=4,
                    num_wells=1,
                    delay=2,
                    reagent_reservoir_volume=1600,#20 * NUM_SAMPLES * 1.1,
                    h_cono=1.95,
                    v_fondo=695,
                    tip_recycling = ['2'])  # Prismatic


    Beads.vol_well = Beads.vol_well_original
    IC.vol_well = IC.vol_well_original
    WashBuffer1.vol_well = WashBuffer1.vol_well_original
    WashBuffer2.vol_well = WashBuffer2.vol_well_original
    ElutionBuffer.vol_well = ElutionBuffer.vol_well_original
    Lysis.vol_well = Lysis.vol_well_original

    ##################
    # Custom functions
    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip,
                       post_dispense=False, post_dispense_vol=20,
                       post_airgap=True, post_airgap_vol=10):
        '''
        x_offset: list with two values. x_offset in source [0] and x_offset in destination [1] i.e. [-1,1]
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
    source_height = 3, post_airgap=True, post_airgap_vol=10,
    post_dispense=False, post_dispense_vol=20,x_offset = [0,0]):
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
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate_mix)
        for _ in range(rounds):
            pipet.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate_mix)
            pipet.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense_mix)
        pipet.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense_mix)
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
        'nest_12_reservoir_15ml', '5', 'Reservoir 12 channel, column 1')

    # Lysis buffer 100ul Deepwell plate and later will be the plate with samples
    ############################################
    magdeck = ctx.load_module('magdeck', '6')
    deepwell = magdeck.load_labware(
        'nest_96_wellplate_2000ul', 'Deepwell plate')


    # Waste plates
    ############################################
    waste_reservoir = ctx.load_labware(
        'nest_1_reservoir_195ml', '9', 'waste reservoir')
    waste = waste_reservoir.wells()[0]  # referenced as reservoir

    # Elution Deepwell plate
    ############################################
    tempdeck = ctx.load_module('tempdeck', '3')
    #tempdeck.set_temperature(temperature)
    qpcr_plate = tempdeck.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt',
        'chilled qPCR final plate')


####################################
    # Load tip_racks
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
               for slot in ['8','10','11']]

    tips300b = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
               for slot in Beads.tip_recycling]
    tips300w1 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
               for slot in WashBuffer1.tip_recycling]
    tips300w2 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
               for slot in WashBuffer2.tip_recycling]

    tips20 = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot, '20µl filter tiprack')
               for slot in ['1']]

################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    WashBuffer1.reagent_reservoir = reagent_res.rows(
    )[0][:1] # position 1

    WashBuffer2.reagent_reservoir = reagent_res.rows(
    )[0][1:2] #position 2

    ElutionBuffer.reagent_reservoir = reagent_res.rows(
    )[0][9:10] #position 10

    Lysis.reagent_reservoir = reagent_res.rows(
    )[0][2:3] #position 3

    IC.reagent_reservoir = reagent_res.rows(
    )[0][10:11] #position 4

    Beads.reagent_reservoir = reagent_res.rows(
    )[0][11:] #position 4

    # columns in destination plates to be filled depending the number of samples
    wb1plate1_destination = deepwell.rows()[0][:num_cols]
    wb2plate1_destination = deepwell.rows()[0][:num_cols]
    qpcr_plate_dest = qpcr_plate.rows()[0][:num_cols]
    sample_plate = deepwell.rows()[0][:num_cols]

    # pipette
    m300 = ctx.load_instrument(
        'p300_multi_gen2', 'right', tip_racks=tips300)  # Load multi pipette

    # pipettes. P1000 currently deactivated
    m20 = ctx.load_instrument(
        'p20_multi_gen2', 'left', tip_racks=tips20)  # Load p300 multi pipette

    # used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0, m20: 0},
        'maxes': {m300: len(tips300)*96, m20: len(tips20)*96}
    }

    ############################################################################
    # STEP 1 Filling with Lysis
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('###############################################')
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        lysis_vol = [100]
        rinse = False  # Only first time

        ########
        # Wash buffer dispense
        if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300b[0].rows()[0][0])
        for i in range(num_cols):
            for j, transfer_vol in enumerate(lysis_vol):
                if (i == 0 and j == 0):
                    rinse = True #Rinse only first transfer
                else:
                    rinse = False
                move_vol_multichannel(m300, reagent = Lysis, source = Lysis.reagent_reservoir[Lysis.col],
                               dest = sample_plate[i], vol = transfer_vol,
                               air_gap_vol = air_gap_vol, x_offset = x_offset,
                               pickup_height = 0.3, rinse = rinse, disp_height = -2,
                               blow_out = True, touch_tip = False, post_airgap=True)

        m300.return_tip()
        tip_track['counts'][m300] += 8
        ctx.comment('Remove Lysis buffer from plate 2')
        ctx.pause('Take plate in 2 to station A and click continue')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 2: TRANSFER IC
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        # Transfer parameters
        start = datetime.now()
        ctx.comment('###############################################')
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        IC_transfer_vol = [ic_vol]
        for i in range(num_cols):
            if not m20.hw_pipette['has_tip']:
                pick_up(m20)
            for j, transfer_vol in enumerate(IC_transfer_vol):
                move_vol_multichannel(m20, reagent=IC, source=IC.reagent_reservoir[IC.col],
                                      dest=sample_plate[i], vol=transfer_vol,
                                      air_gap_vol=air_gap_ic, x_offset=x_offset,
                                      pickup_height=0.2, disp_height = -37.7,
                                      rinse=False, blow_out = True, touch_tip=False, post_airgap=True)

                m20.drop_tip(home_after=False)
                tip_track['counts'][m20] += 8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: PREMIX BEADS
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:

        start = datetime.now()
        ctx.comment('###############################################')
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300b[0].rows()[0][0])
        ctx.comment('Mixing ' + Beads.name)

        # Mixing
        custom_mix(m300, Beads, Beads.reagent_reservoir[Beads.col], vol=120,
                   rounds=10, blow_out=True, mix_height=10, post_dispense=True, source_height=0.3)
        ctx.comment('Finished premixing!')
        ctx.comment('Now, reagents will be transferred to deepwell plate.')
        m300.return_tip()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 4: TRANSFER BEADS and MIX, recycle the tips
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        # Transfer parameters
        start = datetime.now()
        ctx.comment('###############################################')
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        beads_transfer_vol = [beads_vol]  # One round of 20
        rinse = True
        for i in range(num_cols):
            if not m300.hw_pipette['has_tip']:
                m300.pick_up_tip(tips300b[0].rows()[0][i])
            for j, transfer_vol in enumerate(beads_transfer_vol):
                if (i == 0 and j == 0):
                    rinse = True
                else:
                    rinse = False
                # Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(
                    reagent = Beads, cross_section_area = multi_well_rack_area,
                    aspirate_volume = transfer_vol * 8, min_height=0.3, extra_volume=10)
                if change_col == True:  # If we switch column because there is not enough volume left in current reservoir column we mix new column
                    ctx.comment(
                        'Mixing new reservoir column: ' + str(Beads.col))
                    custom_mix(m300, Beads, Beads.reagent_reservoir[Beads.col],
                               vol=120, rounds=10, blow_out=True, mix_height=1,
                               post_dispense=True)
                ctx.comment(
                    'Aspirate from reservoir column: ' + str(Beads.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                move_vol_multichannel(m300, reagent=Beads, source=Beads.reagent_reservoir[Beads.col],
                                      dest=sample_plate[i], vol=transfer_vol,
                                      air_gap_vol=air_gap_vol, x_offset=x_offset,
                                      pickup_height=pickup_height, disp_height = -8,
                                      rinse=rinse, blow_out = True, touch_tip=False, post_airgap=True)

                custom_mix(m300, Beads, sample_plate[i] ,
                                   vol=120, rounds=10, blow_out=True, mix_height=15,
                                   x_offset = x_offset, source_height=0.8, post_dispense=True)
                m300.return_tip(tips300b[0].rows()[0][i])

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 5: INCUBATE WITHOUT MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        # incubate off and on magnet
        magdeck.disengage()
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating OFF magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=2
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 6: INCUBATE WITH MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        magdeck.engage(height=mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=15

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 7: REMOVE SUPERNATANT
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        '''if (wash_buffer2_vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
        # calculate what volume should be transferred in each step
            vol_list=divide_volume(wash_buffer2_vol, pipette_allowed_capacity)'''
        # remove supernatant -> height calculation can be omitted and referred to bottom!
        supernatant_vol = [180]

        for i in range(num_cols):
            offset = find_side(i) * x_offset_rs
            #if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300b[0].rows()[0][i])
            for transfer_vol in supernatant_vol:
                # Pickup_height is fixed here
                pickup_height = 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i + 1))
                ctx.comment('Pickup height is ' +
                            str(pickup_height) + ' (fixed)')
                move_vol_multichannel(m300, reagent=Beads, source=sample_plate[i],
                               dest=waste, vol=transfer_vol, air_gap_vol=air_gap_vol,
                               x_offset=[offset,0], pickup_height=pickup_height, rinse=False,
                               disp_height=-2, blow_out=True, touch_tip=False, post_dispense=True)

            m300.drop_tip(home_after=True)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 8: Washing 1 with Wash1 and MIX
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        magdeck.disengage()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        wash1_vol = [100]
        rinse = WashBuffer1.rinse  # Only first time

        ########
        # W1 isoprop washes
        for i in range(num_cols):
            #if not m300.hw_pipette['has_tip']:
            offset = find_side(i) * - 2 # dispense and mix just over the bead pellet
            m300.pick_up_tip(tips300w1[0].rows()[0][i])
            for j, transfer_vol in enumerate(wash1_vol):
                # Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(
                    WashBuffer1, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from Reservoir column: ' +
                            str(WashBuffer1.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                if i != 0 and j!= 0:
                    rinse = False
                move_vol_multichannel(m300, reagent=WashBuffer1, source=WashBuffer1.reagent_reservoir[WashBuffer1.col],
                               dest=sample_plate[i], vol=transfer_vol, air_gap_vol=air_gap_vol, x_offset=[0,offset],
                               pickup_height=pickup_height, rinse=rinse,
                               disp_height=-4, blow_out=True, touch_tip=False)

                #m300.drop_tip(home_after=True)
                #m300.touch_tip(speed = 20, v_offset = -5)
            custom_mix(m300, reagent=WashBuffer1, location=sample_plate[i], vol=transfer_vol,
                       rounds=2, blow_out=True, mix_height=15, x_offset=[0,offset])
            m300.touch_tip(speed = 20, v_offset = -5)
            m300.return_tip(tips300w1[0].rows()[0][i])

            #tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 9: INCUBATE WITH MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        magdeck.engage(height=mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=15

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 10: REMOVE SUPERNATANT
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        # remove supernatant -> height calculation can be omitted and referred to bottom!
        supernatant_vol = [120]

        for i in range(num_cols):
            offset = find_side(i) * x_offset_rs
            #if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300w1[0].rows()[0][i])
            for transfer_vol in supernatant_vol:
                # Pickup_height is fixed here
                pickup_height = 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i + 1))
                ctx.comment('Pickup height is ' +
                            str(pickup_height) + ' (fixed)')
                move_vol_multichannel(m300, reagent=WashBuffer1, source=sample_plate[i],
                               dest=waste, vol=transfer_vol, air_gap_vol=air_gap_vol,
                               x_offset=[offset,0], pickup_height=pickup_height, rinse=False,
                               disp_height=-2, blow_out=True, touch_tip=True)

            m300.drop_tip(home_after=True)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 11: Washing 2 with Wash2 and MIX
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        magdeck.disengage()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        wash2_vol = [100]
        rinse = WashBuffer2.rinse  # Only first time

        ########
        # isoprop washes
        for i in range(num_cols):
            #if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300w2[0].rows()[0][i])
            offset = find_side(i) * - 2 # dispense and mix just over the bead pellet
            for j, transfer_vol in enumerate(wash2_vol):
                # Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(
                    WashBuffer2, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from Reservoir column: ' +
                            str(WashBuffer1.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                if i != 0 and j!= 0:
                    rinse = False
                move_vol_multichannel(m300, reagent=WashBuffer2, source=WashBuffer2.reagent_reservoir[WashBuffer2.col],
                               dest=sample_plate[i], vol=transfer_vol, air_gap_vol=air_gap_vol, x_offset=[0,offset],
                               pickup_height=pickup_height, rinse=rinse,
                               disp_height=-4, blow_out=True, touch_tip=False)

                #m300.drop_tip(home_after=True)
                m300.touch_tip(speed = 20, v_offset = -5)
            custom_mix(m300, reagent=WashBuffer2, location=sample_plate[i], vol=transfer_vol,
                       rounds=2, blow_out=True, mix_height=15, x_offset=[0,offset])
            m300.touch_tip(speed = 20, v_offset = -5)
            m300.return_tip(tips300w2[0].rows()[0][i])

            #tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 12: INCUBATE WITH MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        magdeck.engage(height=mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=15

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 13: REMOVE SUPERNATANT
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        # remove supernatant -> height calculation can be omitted and referred to bottom!
        supernatant_vol = [120]

        for i in range(num_cols):
            offset = find_side(i) * x_offset_rs
            #if not m300.hw_pipette['has_tip']:
            m300.pick_up_tip(tips300w2[0].rows()[0][i])
            for transfer_vol in supernatant_vol:
                # Pickup_height is fixed here
                pickup_height = 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i + 1))
                ctx.comment('Pickup height is ' +
                            str(pickup_height) + ' (fixed)')
                move_vol_multichannel(m300, reagent=WashBuffer2, source=sample_plate[i],
                               dest=waste, vol=transfer_vol, air_gap_vol=air_gap_vol,
                               x_offset=[offset,0], pickup_height=pickup_height, rinse=False,
                               disp_height=-2, blow_out=True, touch_tip=True)

            m300.drop_tip(home_after=True)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 14: Allow to dry WITH MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        magdeck.engage(height=mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=15

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 15: Transfer Elution buffer ( magnet == OFF)
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        # Water elution
        elut_vol = [elution_vol]
        air_gap_vol_water = 10
        x_offset_w = -1.1
        magdeck.disengage()
        ########
        # Water or elution buffer
        for i in range(num_cols):
            offset = find_side(i)*x_offset_w
            ctx.comment('Side is : '+ str(offset))
            #if not m300.hw_pipette['has_tip']:
            pick_up(m300)
            for transfer_vol in elut_vol:
                # Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(
                    ElutionBuffer, multi_well_rack_area, transfer_vol * 8)
                ctx.comment('Aspirate from Reservoir column: ' + str(ElutionBuffer.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                move_vol_multichannel(m300, reagent=ElutionBuffer, source=ElutionBuffer.reagent_reservoir[ElutionBuffer.col],
                               dest=sample_plate[i], vol=transfer_vol, air_gap_vol=air_gap_vol_water, x_offset=[0,offset],
                               pickup_height=pickup_height, rinse=False, disp_height=-25,
                               blow_out=True, touch_tip=False)

            ctx.comment('Mixing sample with Water and LTA')
            # Mixing
            custom_mix(m300, ElutionBuffer, sample_plate[i], vol=60, rounds=10,
                       blow_out=True, mix_height=155,x_offset=[0,offset],source_height=0.4)
            m300.drop_tip(home_after=True)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 16: INCUBATE WITHOUT MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        # incubate off and on magnet
        magdeck.disengage()
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating OFF magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=2
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 17: INCUBATE WITH MAGNET
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')
        magdeck.engage(height=mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' +
                  format(STEPS[STEP]['wait_time']) + ' seconds.')  # minutes=15

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 19 TRANSFER TO ELUTION PLATE
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()

        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'])
        ctx.comment('###############################################')

        elut_vol_corrected = [x+20 for x in elut_vol] # take 20µl more than needed elution volume
        air_gap_vol_water = 0 # as we will take more volume than wanted an air_gap will be already generated

        for i in range(num_cols):
            offset = find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in elut_vol_corrected:
                # Pickup_height is fixed here
                pickup_height = 0.2
                ctx.comment('Aspirate from deep well column: ' + str(i + 1))
                ctx.comment('Pickup height is ' +
                            str(pickup_height) + ' (fixed)')
                move_vol_multichannel(m300, reagent=Elution, source=sample_plate[i],
                               dest=qpcr_plate_dest[i], vol=transfer_vol, air_gap_vol=air_gap_vol, x_offset=[offset,1],
                               pickup_height=pickup_height, rinse=False,
                               disp_height=-4, blow_out=True, touch_tip=False)

            m300.drop_tip(home_after=True)
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
            #ctx._hw_manager.hardware.set_button_light(1,0,0)
            time.sleep(0.3)
            ctx._hw_manager.hardware.set_lights(rails=True)
            #ctx._hw_manager.hardware.set_button_light(0,0,1)
            time.sleep(0.3)
            ctx._hw_manager.hardware.set_lights(rails=False)
        #ctx._hw_manager.hardware.set_button_light(0,1,0)

    ctx.comment(
        'Finished! \nMove deepwell plates to station C.')
    ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
    ctx.comment('Used racks in total: ' + str(tip_track['counts'][m300] / 96))
