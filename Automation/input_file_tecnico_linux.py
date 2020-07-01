# This file will aim to update and customize the protocol for each sample
# run. Set the number of samples, date, register technician name and create
# the directories to run
#!usr/bin/local/python
# coding=utf-8

from datetime import datetime
import os
import os.path
import pandas as pd
import string
import math
import re
homedir=os.path.expanduser("~")
from multi_well_viral import generate_multi_well_viral
from multi_well_pathogen_IC import generate_multi_well_pathogen_IC
from multi_well_pathogen_R import generate_multi_well_pathogen_R
from multi_mini_well import generate_multi_mini_well

demo_mode=True

# recipes for protocol types [obj. volume per well, allowable remaining nonusable volume in channel]
viral_recipe={'Beads':[20,3],
'Wone':[100,600],
'Wtwo':[100,600],
'IC':[10,3],
'Elution':[50,900],
'Lysis':[100,600],
'MMIX':[20,30]
}

pathogen_recipe={'Beads':[260,600],
'Wone':[300,600],
'Wtwo':[450,600],
'IC':[10,3],
'Elution':[90,600],
'Lysis':[260,600],
'MMIX':[20,30]
}

recipes={'V': viral_recipe, 'P': pathogen_recipe}

main_path = '/home/laboratorio/Documentos/'
code_path = main_path + 'covid19huc/Automation/base_scripts/'
KFV_path = code_path + 'Viral_KF/'
KFP_path = code_path + 'Pathogen_KF/'
excel = main_path + 'covid19huc/Automation/base_scripts/Reference_template.xlsx'

# Function to distinguish between OT and KF protocols
def select_protocol_type(p1, p2):
    ff=False
    while ff==False:
        protocol=input('Introducir protocolo: \nKingfisher VIRAL (V) o Kingfisher PATHOGEN (P) \nProtocolo: ')
        if protocol=='V':
            path=p1 # path
            ff=True
        elif protocol=='P':
            path=p2 # path
            ff=True
        else:
            print('Please, try again')
    return protocol,path

def generate_recipe(mode,cn_samp,recipes,num_samples):
    vol_max_pocillo=12400
    final_recipe={}
    for key in recipes[mode].keys():
        if key not in ['MMIX','Beads','IC']:
            vol_total=math.ceil((recipes[mode][key][0]*cn_samp)/100)*100
            num_cells=math.ceil(recipes[mode][key][0]*cn_samp/vol_max_pocillo)
            vol_pocillo=math.ceil((vol_total/num_cells+recipes[mode][key][1])/100)*100
            final_recipe.update({key: [vol_pocillo,num_cells]})
        elif key == 'MMIX':
            vol_pocillo=recipes[mode][key][0]*(num_samples+2+3)+recipes[mode][key][1]
            num_cells=1
            final_recipe.update({key: [vol_pocillo,num_cells]})
        elif key in ['IC']:
            vol_total=math.ceil((recipes[mode][key][0]*cn_samp)/5)*5
            num_cells=1
            vol_pocillo=math.ceil(((vol_total/num_cells+recipes[mode][key][1])/8/5))*5
            final_recipe.update({key: [vol_pocillo,num_cells]})
        elif (key == 'Beads' and mode == 'V'):
            vol_total=math.ceil((recipes[mode][key][0]*cn_samp)/100)*100
            num_cells=2
            vol_pocillo=math.ceil((vol_total/num_cells+recipes[mode][key][1])/8/10)*10
            final_recipe.update({key: [vol_pocillo,num_cells]})
        elif (key == 'Beads' and mode == 'P'):
            vol_total=math.ceil((recipes[mode][key][0]*cn_samp)/100)*100
            num_cells=math.ceil(recipes[mode][key][0]*cn_samp/vol_max_pocillo)
            vol_pocillo=math.ceil((vol_total/num_cells+recipes[mode][key][1])/100)*100
            final_recipe.update({key: [vol_pocillo,num_cells]})

    return final_recipe

def update_files(final_path,filename,final_data,operation_data):
    f=open(final_path+'/scripts/'+filename, "rt") # open file
    data = f.read()
    for key in final_data.keys():
        if (key in data and key != 'MMIX'):
            data=data.replace('$'+key+'_total_volume',str(final_data[key][0]*final_data[key][1]))
            data=data.replace('$'+key+'_wells',str(final_data[key][1]))
        elif (key in data and key == 'MMIX'):
            data=data.replace('$'+key,str(final_data[key][0]*final_data[key][1]))


    for key in operation_data.keys():
        if key in data:
            data=data.replace(key,str(operation_data[key]))

    f=open(os.path.join(final_path+'/scripts/',filename), "wt")
    f.write(data)
    f.close()

def update_readme(final_path,filename,protocol,imagepath,operation_data):
    f=open(os.path.join(final_path,filename), "rt") # open file
    data = f.read()
    for key in operation_data.keys():
        if key in data:
            data=data.replace(key,str(operation_data[key]))

    if protocol == 'V':
        data=data.replace('img_src',str(imagepath[0]))
        data=data.replace('img_mmwell',str(imagepath[1]))


    if protocol == 'P':
        data=data.replace('img_src_R',str(imagepath[0]))
        data=data.replace('img_src_IC',str(imagepath[1]))
        data=data.replace('img_mmwell',str(imagepath[2]))

    f=open(os.path.join(final_path,filename), "wt")
    f.write(data)
    f.close()

###############################################################################
def main():
    if demo_mode==False:
        # Read the excel file from the run and obtain the dictionary of samples
        excel = '/home/laboratorio/Escritorio/prueba.xlsx'
        df = pd.read_excel (excel,
         sheet_name='Deepwell layout', header = None, index_col = 0)
        df = df.iloc[2:]
        df_dict = df.to_dict('index')
        merged_dict={}
        for key in df_dict:
            for key2 in df_dict[key]:
                merged_dict[str(key)+format(key2)]=df_dict[key][key2]

        # count number of declared elements in Dictionary
        num_samples_control = 0
        for elem in merged_dict.values():
            if elem != 0:
                num_samples_control += 1

        # Get sample data from user
        control=False
        while control==False:
            num_samples = int(input('Número de muestras a procesar (excluidos PC + NC): '))
            if (num_samples>0 and num_samples<=94):
                control=True
            else:
                print('Número de muestras debe ser un número entre 1 y 94 (2 son los controles)')
                print('------------------------------------------------------------------------')
        print('El número de muestras registradas en el excel es: '+str(num_samples_control))
        if num_samples_control!=num_samples:
            print('Error: El número de muestras entre excel y reportado no coincide, revisar por favor.')
            exit()
        else:
            print('El número de muestras coincide')
            print('------------------------------------------------------------------------')

        # Get technician name
        control=False
        while control==False:
            tec_name = '\''+(input('Nombre del técnico (usuario): '))+'\''
            print('------------------------------------------------------------------------')
            if isinstance(tec_name, str):
                control=True
            else:
                print('Introduce tu usuario HUC, por favor')

        # Get run session ID
        control=False
        while control==False:
            id = int(input('ID run: '))
            print('------------------------------------------------------------------------')
            if isinstance(id,int):
                control=True
            else:
                print('Por favor, assigna un ID numérico para este RUN')
    else:
        # Get sample data from user
        control=False
        while control==False:
            num_samples = int(input('Número de muestras a procesar (excluidos PC + NC): '))
            if (num_samples>0 and num_samples<=94):
                control=True
        id=43001
        excel = '/home/laboratorio/Escritorio/prueba.xlsx'
        tec_name='demo'


    # Get date
    fecha=datetime.now()
    t_registro='\''+fecha.strftime("%m/%d/%Y, %H:%M:%S")+'\''
    h_registro=fecha.strftime("%H:%M")
    dia_registro=fecha.strftime("%Y_%m_%d")

    # select the type of protocol to be run
    [protocol,protocol_path]=select_protocol_type(KFV_path, KFP_path)
    print('------------------------------------------------------------------------')

    num_samples_c = math.ceil(num_samples/8)*8 # corrected num_samples value to calculate needed volumes
    num_cols = math.ceil(num_samples_c/8)
    final_data=generate_recipe(protocol,num_samples_c,recipes,num_samples)
    print(final_data)
    operation_data={'$technician': str(tec_name), '$num_samples': str(num_samples),
                '$date': str(t_registro), '$run_id': str(id),
                '$hora': str(h_registro), '$dia': str(dia_registro),
                '$num_s_corrected': str(num_samples_c),
                '$num_cols': str(num_cols),
                '$MMIX': str(final_data['MMIX'][0])
                }

    #determine output path
    run_name = str(dia_registro)+'_OT'+str(id)+'_'+protocol
    final_path=os.path.join(main_path,run_name)

    # create folder directory in case it doesn't already exist and copy excel registry file there

    if not os.path.isdir(final_path):
        os.mkdir(final_path)
        os.mkdir(final_path+'/scripts')
        os.mkdir(final_path+'/results')
        os.mkdir(final_path+'/logs')
        os.system('cp ' + excel +' '+ final_path+'/OT'+str(id)+'_samples.xlsx') # copy excel input file to destination
        if protocol == 'V':
            os.system('cp ' +main_path +'covid19huc/Automation/volumes_viral_readme.html' + ' ' + final_path + '/readme.html')
            pV=generate_multi_well_viral(final_path,final_data)
            mini_well=generate_multi_mini_well(final_path,final_data,protocol)
            update_readme(final_path,'readme.html',protocol,[pV,mini_well],operation_data)
        elif protocol == 'P':
            os.system('cp ' +main_path +'covid19huc/Automation/volumes_pathogen_readme.html' + ' ' + final_path + '/readme.html')
            pB=generate_multi_well_pathogen_IC(final_path,final_data)
            pR=generate_multi_well_pathogen_R(final_path,final_data)
            mini_well=generate_multi_mini_well(final_path,final_data,protocol)
            update_readme(final_path,'readme.html',protocol,[pR,pB,mini_well],operation_data)

    else:
        print('BEWARE! This protocol and ID run already exists! Exitting...')
        exit()


    # move protocol .py files to final destination
    for file in os.listdir(protocol_path): # look for all protocols in folder
        if file.endswith('.py'):
            position=re.search('_template',file).start() # find _ position after the name and get value
            filename=file[:position]+'_'+str(dia_registro)+'_OT'+str(id)+'.py' # assign a filename date + station name + id
            os.system('cp ' + os.path.join(protocol_path,file) +' '+ os.path.join(final_path+'/scripts/',filename))

    # change values to protocols for final user
    for filename in os.listdir(final_path+'/scripts/'):
        if re.match('KA',filename):
            update_files(final_path,filename,final_data,operation_data)

        elif re.match('KB',filename):
            update_files(final_path,filename,final_data,operation_data)

        elif re.match('KC',filename):
            update_files(final_path,filename,final_data,operation_data)

        else:
            print('No files found')

if __name__ == '__main__':
    main()

    print('Success!')
