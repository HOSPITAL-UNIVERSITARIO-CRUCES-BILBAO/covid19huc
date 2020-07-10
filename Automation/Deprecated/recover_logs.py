import os
import time

run_id ='4300'
run_type ='pathogen'

logs_path='/home/downloads'

rpi_path = os.path.join('/var/lib/jupyter/notebooks/',run_id)
main_path = os.path.join('/home/downloads',run_id)

# filenames as they are in the robot log
file_stationa = '/KA_SampleSetup_'+run_type+'_time_log.txt'

file_stationb_beads = '/KB_beads_'+run_type+'_time_log.txt'

file_stationb_plates = '/KB_PlateFilling_'+run_type+'_time_log.txt'

file_stationc = '/KC_qPCR_'+run_type+'_time_log.txt'

filenames = [file_stationa,file_stationb_beads,file_stationb_plates,file_stationc]


# IP routing for these robots
IP_list = ['10.39.219.'+str(i) for i in [12,14,16]]

if not os.path.isdir(main_path):
    os.mkdir(main_path)
    os.mkdir(main_path+'/logs')

for IP,filename in zip(IP_list,filenames):
    try:
        os.system('scp root@'+IP+':'+os.path.join(rpi_path,filename)+' '+str(os.path.join(main_path,filename)))
        time.sleep(2)
    except:
        print('file '+filename+' could not be exported')
