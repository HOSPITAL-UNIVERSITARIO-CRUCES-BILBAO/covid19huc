
def generate_multi_well_pathogen_IC(path, recipe):
    import matplotlib.pyplot as plt
    import math

    ###################
    n_beads_max = 3  # max number of columns with beads
    n_lysis_max = 3  # max number of columns with lysis
    color_one = 'yellowgreen'
    color_two = 'indianred'
    bar_label = []
    ###################

    plt.rcParams.update({'font.size': 16})
    #plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
    x = [i + 1 for i in range(12)]
    val = 10
    y = [val for i in range(12)]
    toplabels = ['IC']
    toplabels.extend(['' for i in range(4)])
    toplabels.extend(['', 'Beads', '', '', '', 'Lysis', ''])

    # generate color table
    # first for IC
    colors = [color_one]
    bar_label.extend([recipe['IC'][0] for i in range(recipe['IC'][1])])

    # second for 4 white positions
    colors.extend(['white' for i in range(4)])
    bar_label.extend([0 for i in range(4)])

    # now beads
    if recipe['Beads'][1] == n_beads_max:
        colors.extend([color_two for i in range(recipe['Beads'][1])])
        # generate labels of how much volume should be put into each well
        bar_label.extend([recipe['Beads'][0] for i in range(recipe['Beads'][1])])
    else:
        colors.extend([color_two for i in range(recipe['Beads'][1])])
        colors.extend(['white' for i in range(
            n_beads_max - recipe['Beads'][1])])
        bar_label.extend([recipe['Beads'][0] for i in range(recipe['Beads'][1])])
        bar_label.extend([0 for i in range(n_beads_max - recipe['Beads'][1])])

    # a white one
    colors.extend(['white'])
    bar_label.extend([0])

    # finally lysis
    if recipe['Lysis'][1] == n_lysis_max:
        colors.extend([color_one for i in range(recipe['Lysis'][1])])
        bar_label.extend([recipe['Lysis'][0] for i in range(recipe['Lysis'][1])])
    else:
        colors.extend([color_one for i in range(recipe['Lysis'][1])])
        colors.extend(['white' for i in range(
            n_lysis_max - recipe['Lysis'][1])])
        bar_label.extend([recipe['Lysis'][0]  for i in range(recipe['Lysis'][1])])
        bar_label.extend([0 for i in range(n_lysis_max - recipe['Lysis'][1])])
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
    plt.title('IC Reservoir')

    # plt.show()
    plt.savefig(path + '/pathogen_IC_multi_well_layout.png',
                dpi=300, bbox_inches='tight')
    return path + '/pathogen_IC_multi_well_layout.png'
