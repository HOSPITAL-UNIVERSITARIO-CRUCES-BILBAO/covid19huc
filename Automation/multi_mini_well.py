def generate_multi_mini_well(path, recipe,mode,num_cols):
    import matplotlib.pyplot as plt
    import math
    import numpy as np


    plt.rcParams.update({'font.size': 16})
    #plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
    color_one='yellowgreen'
    color_two='indianred'

    x = [i+1 for i in range(12)]
    val = 10
    y = [val for i in range(12)]


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks([])
    axT=ax.twiny()
    plt.ylim([-0.5,7.5])
    plt.xlim([-1,14])



    #create plot depending on protocol type
    if mode=='V':
        #assign top labels
        toplabels=['IC','Be','Be']
        toplabels.extend(['' for i in range(4,13)])
        #create wells to fill
        barplot=ax.scatter([1 for x in range(8)],np.arange(8),s=200,c = 'g')
        #create empty wells
        for i in range(4,13):
            barplot=ax.scatter([i for x in range(8)],np.arange(8),s=200, facecolors='none', edgecolors='k')
        # colour beads _wells
        if num_cols==1:
            barplot=ax.scatter([2 for x in range(8)],np.arange(8),s=200,c = 'r')
            barplot=ax.scatter([3 for x in range(8)],np.arange(8),s=200,facecolors='none', edgecolors='k')
        else:
            barplot=ax.scatter([2 for x in range(8)],np.arange(8),s=200,c = 'r')
            barplot=ax.scatter([3 for x in range(8)],np.arange(8),s=200,c = 'r')
        #create empty wells
        for i in range(4,13):
            barplot=ax.scatter([i for x in range(8)],np.arange(8),s=200, facecolors='none', edgecolors='k')
        #create annotations
        dot_label = [recipe['IC'][0]]
        if (num_cols % 2 == 0 and num_cols>1):
            dot_label.extend([recipe['Beads'][0] for i in range(recipe['Beads'][1]) ])
        elif num_cols == 1:
            dot_label.extend([recipe['Beads'][0] + 10])
            dot_label.extend([0])
        else:
            dot_label.extend([recipe['Beads'][0] + 10])
            dot_label.extend([recipe['Beads'][0] - 10])
        dot_label = [str(i)+'µl' for i in dot_label]

    elif mode == 'P':
        # assign top labels
        toplabels=['IC']
        toplabels.extend(['' for i in range(2,13)])
        #create wells to fill
        barplot=ax.scatter([1 for x in range(8)],np.arange(8),s=200,c = 'g')
        #create empty wells
        for i in range(2,13):
            barplot=ax.scatter([i for x in range(8)],np.arange(8),s=200, facecolors='none', edgecolors='k')
        #create empty wells
        for i in range(2,13):
            barplot=ax.scatter([i for x in range(8)],np.arange(8),s=200, facecolors='none', edgecolors='k')
        #create annotations
        dot_label = [recipe['IC'][0]]
        dot_label = [str(i)+'µl' for i in dot_label]

    # well annotation parameters
    px=[1,2,3]
    py=[6.8,6.8,6.8]
    for i, txt in enumerate(dot_label):
        ax.annotate(txt, (px[i], py[i]),rotation=30)

    axT.set_xlim(ax.get_xlim())
    ax.set_xticks(x)
    axT.set_xticks(x)
    axT.set_xticklabels(toplabels)
    plt.title('Internal control and beads* reservoir')
    #lt.show()

    plt.savefig(path+'/multi_mini_well_layout.png',dpi=300,bbox_inches='tight')
    return path+'/multi_mini_well_layout.png'
