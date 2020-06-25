def generate_multi_well_pathogen_R(path, recipe):
    import matplotlib.pyplot as plt
    import math

    ###################
    wash_two_max = 4  # max number of columns with beads
    wash_one_max = 3  # max number of columns with lysis
    color_one = 'yellowgreen'
    color_two = 'indianred'
    bar_label = []
    ###################

    plt.rcParams.update({'font.size': 16})
    #plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
    x = [i+1 for i in range(12)]
    val = 10
    y = [val for i in range(12)]
    toplabels=['Elution']
    toplabels.extend(['' for i in range(4)])
    toplabels.extend(['','Wash 2     ','','','','Wash 1',''])

    # generate color table
    # first for IC
    colors = [color_one]
    bar_label.extend([recipe['Elution'][0] for i in range(recipe['Elution'][1])])

    # second for 4 white positions
    colors.extend(['white' for i in range(3)])
    bar_label.extend([0 for i in range(3)])

    # now wash two

    if recipe['Wtwo'][1] == wash_two_max:
        colors.extend([color_two for i in range(recipe['Wtwo'][1])])
        # generate labels of how much volume should be put into each well
        bar_label.extend([recipe['Wtwo'][0] for i in range(recipe['Wtwo'][1])])
    else:
        colors.extend([color_two for i in range(recipe['Wtwo'][1])])
        colors.extend(['white' for i in range(wash_two_max-recipe['Wtwo'][1])])
        bar_label.extend([recipe['Wtwo'][0] for i in range(recipe['Wtwo'][1])])
        bar_label.extend([0 for i in range(wash_two_max - recipe['Wtwo'][1])])

    # a white one
    colors.extend(['white'])
    bar_label.extend([0])

    # finally wash one
    if recipe['Wone'][1] == wash_two_max:
        colors.extend([color_one for i in range(recipe['Wone'][1])])
        # generate labels of how much volume should be put into each well
        bar_label.extend([recipe['Wone'][0] for i in range(recipe['Wone'][1])])
    else:
        colors.extend([color_one for i in range(recipe['Wone'][1])])
        colors.extend(['white' for i in range(wash_one_max-recipe['Wone'][1])])
        bar_label.extend([recipe['Wone'][0] for i in range(recipe['Wone'][1])])
        bar_label.extend([0 for i in range(wash_one_max - recipe['Wone'][1])])
    ########
    bar_label = [str(i) + ' µl' for i in bar_label]

    #######
    #generate plot frame
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks([])
    axT = ax.twiny()

    # generate plot
    barplot = ax.bar(x, y)
    plt.ylim([-0.5, 10.5])

        # add volumes to each well for each bar
    for idx, rect in enumerate(barplot):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 0.5 * height,
                bar_label[idx],
                ha='center', va='center', rotation=90)

    for i in range(len(barplot)):
        barplot[i].set_color(colors[i])
        if colors[i] == 'white':
            barplot[i].set_ec('k')

    axT.set_xlim(ax.get_xlim())
    ax.set_xticks(x)
    axT.set_xticks(x)
    axT.set_xticklabels(toplabels)
    plt.title('R Reservoir')

    # plt.show()
    plt.savefig(path + '/pathogen_R_multi_well_layout.png',
                dpi=300, bbox_inches='tight')

'''import matplotlib.pyplot as plt

###################
wash_two_max = 4  # max number of columns with beads
wash_one_max = 3  # max number of columns with lysis
color_one = 'yellowgreen'
color_two = 'indianred'
bar_label = []
###################

plt.rcParams.update({'font.size': 16})
#plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
x = [i+1 for i in range(12)]
val = 10
y = [val for i in range(12)]
toplabels=['Elution']
toplabels.extend(['' for i in range(4)])
toplabels.extend(['','Wash 2     ','','','','Wash 1',''])

#### generate color table
colors=[color_one]
colors.extend(['white' for i in range(3)])




bar_label=[1000,0,0,0,2000,0,0,0,0,1000,0,0]

bar_label = [str(i)+'µl' for i in bar_label]

if wash_one == 3:
    colors.extend([color_one for i in range(wash_one)])
else:
    colors.extend([color_one for i in range(wash_one)])
    colors.extend(['white' for i in range(3-wash_one)])
########

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks([])
axT=ax.twiny()

barplot=ax.bar(x,y)
plt.ylim([-0.5,10.5])

# set color to bar
for i in range(len(barplot)):
    barplot[i].set_color(colors[i])
    if colors[i]=='white':
        barplot[i].set_ec('k')

# add volumes to each well for each bar
for idx,rect in enumerate(barplot):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height,
                bar_label[idx],
                ha='center', va='center', rotation=90)

axT.set_xlim(ax.get_xlim())
ax.set_xticks(x)
axT.set_xticks(x)
axT.set_xticklabels(toplabels)

plt.show()
'''
