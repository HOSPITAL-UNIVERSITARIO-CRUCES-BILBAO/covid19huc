
def generate_multi_well(path):
    import os
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 16})
    #plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
    color_one='yellowgreen'
    color_two='indianred'

    x = [i+1 for i in range(12)]
    val = 10
    y = [val for i in range(12)]
    toplabels=['W1','W2','L']
    toplabels.extend(['' for i in range(3,9)])
    toplabels.extend(['E','IC','B'])

    bar_label=[1000,0,0,0,2000,0,0,0,0,1000,0,0]
    bar_label = [str(i)+'µl' for i in bar_label]

    colors=[color_one, color_two,color_one]
    colors.extend(['white' for i in range(3,9)])
    colors.extend([color_one, color_two,color_one])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks([])
    axT=ax.twiny()

    barplot=ax.bar(x,y)
    plt.ylim([-0.5,10.5])

    for i in range(len(barplot)):
        barplot[i].set_color(colors[i])

    for i in range(3,9):
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
    plt.title('Reservoir')
    #plt.show()
    plt.savefig(path+'/viral_multi_well_layout.png',dpi=300,bbox_inches='tight')
