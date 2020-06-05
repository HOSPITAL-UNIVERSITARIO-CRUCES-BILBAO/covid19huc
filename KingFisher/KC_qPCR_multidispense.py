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
    'protocolName': 'Station C Kingfisher Pathogen qPCR setup Version 2',
    'author': 'Malen Aguirregabiria, Aitor Gastaminza & José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Clínic Barcelona, Hospital Universitario Cruces Bilbao',
    'apiLevel': '2.0',
    'description': 'Protocol for Kingfisher sample setup (C) - Viral Kit (ref )'

}
'''
'technician': '$technician',
'date': '$date'
'''
#Defined variables
##################
NUM_SAMPLES = 16
#NUM_SAMPLES = NUM_SAMPLES -1 #Remove last sample (PC), done manually


air_gap_vol = 30
air_gap_sample = 2
run_id = '$run_id'

# Tune variables
size_transfer = 4  # Number of wells the distribute function will fill
volume_mmix = 30  # Volume of transfered master mix
volume_sample = 5  # Volume of the sample
volume_mmix_available = 1000 #(NUM_SAMPLES * 1.5 * volume_mmix)  # Total volume of first screwcap
extra_dispensal = 10  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.25  # Diameter of the screwcap
temperature = 25  # Temperature of temp module
volume_cone = 50  # Volume in ul that fit in the screwcap cone
x_offset = [0,0]
pipette_allowed_capacity=170
#size_transfer = math.floor(pipette_allowed_capacity / volume_mmix)

#############################################################
#Available master mastermixes
#############################################################
MMIX_available={1: 'Seegene', 2: 'Universal', 3: 'Universal_IDT',4: 'Clinic', 5: 'Multiplex'}

mmix_selection = 5 # select the mastermix to be used

MMIX_vol={1: [17,1], 2: [20,1], 3: [20,1], 4: [40,2], 5:[20,2]} # volume of mastermixes per sample and number of wells in which is distributed
MMIX_recipe={1: [5, 5, 5, 2], 2: [8, 5, 1, 2, 2, 1, 1], 3: [12, 5, 1, 1, 1], 4: [1], 5:[6.25,1.25,12.5]} # Reactive volumes for the mmix

size_transfer = math.floor(pipette_allowed_capacity / MMIX_vol[mmix_selection][0]) # Number of wells the distribute function will fill

MMIX_make_location = 9 # Cell C1 in which the first tube for the MMIX will be placed

volume_mmix = MMIX_vol[mmix_selection][0]  # Volume of transfered master mix per well

MMIX_make = {}
for mmix_type in MMIX_recipe.keys():
    MMIX_make[mmix_type] = []
    for needed_vol in MMIX_recipe[mmix_type]:
        MMIX_make[mmix_type].append(needed_vol * NUM_SAMPLES * 1.1)

volume_mmix_available = (NUM_SAMPLES * 1.1 * MMIX_vol[mmix_selection][0])  # Total volume of mastermix that will be prepared


#############################################



# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):
    import os
    #from opentrons.drivers.rpi_drivers import gpio
    #gpio.set_rail_lights(False) #Turn off lights (termosensible reagents)
    ctx.comment('Actual used columns: ' + str(num_cols))

    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': True, 'description': 'Make MMIX'},
        2: {'Execute': False, 'description': 'Transfer MMIX'},
        3: {'Execute': False, 'description': 'Transfer elution'}
    }

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    #Folder and file_path for log time
    folder_path = '/var/lib/jupyter/notebooks/'+run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/KC_qPCR_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none'):
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

    # Reagents and their characteristics
    MMIX = Reagent(name = 'Master Mix',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = volume_mmix_available,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    MMIX_components = Reagent(name = 'MMIX_component',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    Elution = Reagent(name='Elution',
                      rinse=False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume=50,
                      delay=0,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=0,
                      v_fondo=0
                      )

    Tackpath = Reagent(name = 'MMIX_multiplex_tackpath',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    covid_assay = Reagent(name = 'Covid19_assay',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    mmix_water = Reagent(name = 'nuclease_free_water',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    MMIX.vol_well = MMIX.vol_well_original
    Elution.vol_well = Elution.vol_well_original
    Tackpath.vol_well = Tackpath.vol_well_original
    covid_assay.vol_well = covid_assay.vol_well_original
    mmix_water.vol_well = mmix_water.vol_well_original

    MMIX_components=[Tackpath,covid_assay,mmix_water]


    ##################
    # Custom functions
    def divide_volume(volume,max_vol):
        num_transfers=math.ceil(volume/max_vol)
        vol_roundup=math.ceil(volume/num_transfers)
        last_vol = volume - vol_roundup*(num_transfers-1)
        vol_list = [vol_roundup for v in range(1,num_transfers)]
        vol_list.append(last_vol)
        return vol_list

    def divide_destinations(l, n):
        a=[]
        # Divide the list of destinations in size n lists.
        for i in range(0, len(l), n):
            a.append( l[i:i + n])

        return a

    def distribute_custom(pipette, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        air_gap=10
        pipette.aspirate((len(dest) * volume) +
                         extra_dispensal, src.bottom(pickup_height), rate = reagent.flow_rate_aspirate)
        pipette.touch_tip(speed=20, v_offset=-5)
        pipette.move_to(src.top(z=5))
        pipette.aspirate(air_gap, rate = reagent.flow_rate_aspirate)  # air gap

        for d in dest:
            pipette.dispense(air_gap, d.top(), rate = reagent.flow_rate_dispense)
            drop = d.top(z = disp_height)
            pipette.dispense(volume, drop, rate = reagent.flow_rate_dispense)
            ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent
            pipette.move_to(d.top(z=5))
            pipette.aspirate(air_gap,d.top(z=5), rate = reagent.flow_rate_aspirate)  # air gap

        try:
            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pipette.blow_out(waste_pool.bottom(pickup_height + 3))
        return (len(dest) * volume)

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
            pipet.dispense(post_airgap_vol, dest.top(z = 5), rate = 2)



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

    def calc_height(reagent, cross_section_area, aspirate_volume, min_height = 0.5, extra_volume = 50):
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

    ####################################
    # load labware and modules
    # 24 well rack
    tuberack = ctx.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap', '1',
        'Bloque Aluminio opentrons 24 screwcaps 2000 µL ')

    ############################################
    # tempdeck
    tempdeck = ctx.load_module('tempdeck', '4')
    #tempdeck.set_temperature(temperature)

    ##################################
    # qPCR plate - final plate, goes to PCR
    qpcr_plate = tempdeck.load_labware(
        'abi_fast_qpcr_96_alum_opentrons_100ul',
        'chilled qPCR final plate')

    ##################################
    # Sample plate - comes from B
    source_plate = ctx.load_labware(
        "kingfisher_std_96_wellplate_550ul", '2',
        'chilled KF plate with elutions (alum opentrons)')
    samples = source_plate.wells()[:NUM_SAMPLES]

    ##################################
    # Load Tipracks
    tips20 = [
        ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['5']
    ]

    tips200 = [
        ctx.load_labware('opentrons_96_filtertiprack_200ul', slot)
        for slot in ['6']
    ]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    MMIX.reagent_reservoir = tuberack.rows()[0][:MMIX.num_wells] # 1 row, 2 columns (first ones)
    MMIX_components_location=tuberack.wells()[MMIX_make_location:(MMIX_make_location + len(MMIX_make[mmix_selection]))]
    ctx.comment('Wells in: '+ str(tuberack.rows()[0][:MMIX.num_wells]) + ' element: '+str(MMIX.reagent_reservoir[MMIX.col]))

    # setup up sample sources and destinations
    samples = source_plate.wells()[:NUM_SAMPLES]
    samples_multi = source_plate.rows()[0][:num_cols]
    pcr_wells = qpcr_plate.wells()[:(NUM_SAMPLES)]
    pcr_wells_multi = qpcr_plate.rows()[0][:num_cols]
    # Divide destination wells in small groups for P300 pipette
    dests = list(divide_destinations(pcr_wells, size_transfer))


    # pipettes
    m20 = ctx.load_instrument(
        'p20_multi_gen2', mount='left', tip_racks=tips20)
    p300 = ctx.load_instrument(
        'p300_single_gen2', mount='right', tip_racks=tips200)

    # used tips counter
    tip_track = {
        'counts': {p300: 0, m20: 0},
        'maxes': {p300: len(tips200) * 96, m20: len(tips20)*96}
    }

    ########################################################################
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        if not ctx.is_simulating():
            if tip_track['counts'][pip] == tip_track['maxes'][pip]:
                ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0

        if not pip.hw_pipette['has_tip']:
            pip.pick_up_tip()
    #########################################################################

    ############################################################################
    # STEP 1: Make Master MIX
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        # Check if among the pipettes, p300_single is installed
        used_vol=[]
        for i,[source, vol] in enumerate(zip(MMIX_components_location, MMIX_make[mmix_selection])):
            pick_up(p300)
            if (vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
            # calculate what volume should be transferred in each step
                vol_list=divide_volume(vol, pipette_allowed_capacity)
                for vol in vol_list:
                    move_vol_multichannel(p300, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                    vol=vol, air_gap_vol=air_gap_vol, x_offset = x_offset,pickup_height=1,
                    rinse=False, disp_height=-10,blow_out=True, touch_tip=True)
            else:
                move_vol_multichannel(p300, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset,pickup_height=1,
                rinse=False, disp_height=-10,blow_out=True, touch_tip=True)

            if i<len(MMIX_components):
                p300.drop_tip()
            else:
                custom_mix(p300, reagent = MMIX, location = MMIX.reagent_reservoir[0], vol = 180, rounds = 3,
                blow_out = True, mix_height = 2, x_offset = x_offset)
                p300.drop_tip()

            tip_track['counts'][p300]+=1



        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 2: Transfer Master MIX
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=True) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p300.pick_up_tip()

        used_vol=[]
        for dest in dests:
            aspirate_volume=volume_mmix * len(dest) + extra_dispensal
            [pickup_height,col_change]=calc_height(MMIX, area_section_screwcap, aspirate_volume)
            ctx.comment('destinations '+ str(dest))
            used_vol_temp = distribute_custom(p300, reagent= MMIX, volume = volume_mmix,
            src = MMIX.reagent_reservoir[MMIX.col], dest = dest,
            waste_pool = MMIX.reagent_reservoir[MMIX.col], pickup_height = pickup_height,
            extra_dispensal = extra_dispensal)

            used_vol.append(used_vol_temp)

        p300.drop_tip()
        tip_track['counts'][p300]+=1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: TRANSFER Samples
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('pcr_wells')
        #Loop over defined wells
        for s, d in zip(samples_multi, pcr_wells_multi):
            m20.pick_up_tip()
            #Source samples
            move_vol_multichannel(m20, reagent = Elution, source = s, dest = d,
            vol = volume_sample, air_gap_vol = air_gap_sample, x_offset = x_offset,
                   pickup_height = 0.2, disp_height = -10, rinse = False,
                   blow_out=True, touch_tip=True, post_airgap=True)
            m20.drop_tip()
            tip_track['counts'][m20]+=8

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

    time.sleep(2)
    import os
    #os.system('mpg123 -f -8000 /etc/audio/speaker-test.mp3 &')

    if STEPS[1]['Execute'] == True:
        total_used_vol = np.sum(used_vol)
        total_needed_volume = total_used_vol
        ctx.comment('Total Master Mix used volume is: ' + str(total_used_vol) + '\u03BCl.')
        ctx.comment('Needed Master Mix volume is ' +
                    str(total_needed_volume + extra_dispensal*len(dests)) +'\u03BCl')
        ctx.comment('Used Master Mix volumes per run are: ' + str(used_vol) + '\u03BCl.')
        ctx.comment('Master Mix Volume remaining in tubes is: ' +
                    format(np.sum(MMIX.unused)+extra_dispensal*len(dests)+MMIX.vol_well) + '\u03BCl.')
        ctx.comment('200 ul Used tips in total: ' + str(tip_track['counts'][p300]))
        ctx.comment('200 ul Used racks in total: ' + str(tip_track['counts'][p300] / 96))

    if STEPS[2]['Execute'] == True:
        ctx.comment('20 ul Used tips in total: ' + str(tip_track['counts'][m20]))
        ctx.comment('20 ul Used racks in total: ' + str(tip_track['counts'][m20] / 96))

    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        #ctx._hw_manager.hardware.set_button_light(1,0,0)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        #ctx._hw_manager.hardware.set_button_light(0,0,1)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=False)
    #ctx._hw_manager.hardware.set_button_light(0,1,0)

    ctx.comment('Finished! \nMove plate to PCR')
