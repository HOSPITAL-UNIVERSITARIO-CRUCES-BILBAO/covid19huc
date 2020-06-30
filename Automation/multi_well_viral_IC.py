import matplotlib.pyplot as plt
import math
import numpy as np

def generate_recipe(mode,cn_samp,recipes,num_samples):
    vol_max_pocillo=12400
    final_recipe={}
    for key in recipes[mode].keys():
        if key != 'MMIX':
            vol_total=math.ceil((recipes[mode][key][0]*cn_samp)/100)*100
            num_cells=math.ceil(recipes[mode][key][0]*cn_samp/vol_max_pocillo)
            vol_pocillo=math.ceil((vol_total/num_cells+recipes[mode][key][1])/100)*100
            final_recipe.update({key: [vol_pocillo,num_cells]})
        elif key == 'MMIX':
            vol_pocillo=recipes[mode][key][0]*(num_samples+2+3)+recipes[mode][key][1]
            num_cells=1
            final_recipe.update({key: [vol_pocillo,num_cells]})
    return final_recipe

viral_recipe={'Beads':[20,800],
'Wone':[100,600],
'Wtwo':[100,600],
'IC':[10,1500],
'Elution':[50,900],
'Lysis':[100,600],
'MMIX':[20,30]
}
pathogen_recipe={'Beads':[260,600],
'Wone':[300,600],
'Wtwo':[450,600],
'IC':[10,1500],
'Elution':[90,600],
'Lysis':[260,600],
'MMIX':[20,30]
}

recipes={'V': viral_recipe, 'P': pathogen_recipe}
num_samples=96
num_samples_c = math.ceil(num_samples/8)*8 # corrected num_samples value to calculate needed volumes
num_cols = math.ceil(num_samples_c/12)

recipe=generate_recipe('V',num_samples_c,recipes,num_samples)


plt.rcParams.update({'font.size': 16})
#plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
color_one='yellowgreen'
color_two='indianred'

x = [i+1 for i in range(12)]
val = 10
y = [val for i in range(12)]
toplabels=['IC','B','B']
toplabels.extend(['' for i in range(4,13)])

bar_label=[recipe['Wone'][0],recipe['Wtwo'][0],recipe['Lysis'][0],0,0,0,0,0,0,recipe['Elution'][0],recipe['IC'][0],recipe['Beads'][0]]
bar_label = [str(i)+'µl' for i in bar_label]

colors=[color_one, color_two,color_one]
colors.extend(['white' for i in range(3,9)])
colors.extend([color_one, color_two,color_one])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks([])
axT=ax.twiny()

barplot=ax.scatter([1 for x in range(8)],np.arange(8),s=200,c = 'g')
barplot=ax.scatter([2 for x in range(8)],np.arange(8),s=200,c = 'r')
barplot=ax.scatter([3 for x in range(8)],np.arange(8),s=200,c = 'r')
for i in range(4,13):
    barplot=ax.scatter([i for x in range(8)],np.arange(8),s=200, facecolors='none', edgecolors='k')
plt.ylim([-0.5,7.5])
plt.xlim([-1,14])

px=[1,2,3]
py=[7,6,7]

dot_label = [125, 125, 125]
dot_label = [str(i)+'µl' for i in dot_label]

for i, txt in enumerate(dot_label):
    ax.annotate(txt, (px[i], py[i]))

'''for i in range(len(barplot)):
    barplot[i].set_color(colors[i])

for i in range(3,9):
    barplot[i].set_ec('k')

# add volumes to each well for each bar
for idx,rect in enumerate(barplot):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height,
                bar_label[idx],
                ha='center', va='center', rotation=90)'''

axT.set_xlim(ax.get_xlim())
ax.set_xticks(x)
axT.set_xticks(x)
axT.set_xticklabels(toplabels)
plt.title('Reservoir')
plt.show()
    #plt.savefig(path+'/viral_multi_well_layout.png',dpi=300,bbox_inches='tight')
