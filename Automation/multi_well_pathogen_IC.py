import matplotlib.pyplot as plt

###################
n_beads=3
n_lysis=3
color_one='yellowgreen'
color_two='indianred'
###################

plt.rcParams.update({'font.size': 16})
#plt.rcParams.update({'font.family': 'Computer Modern Sans Serif'})
x = [i+1 for i in range(12)]
val = 10
y = [val for i in range(12)]
toplabels=['IC']
toplabels.extend(['' for i in range(4)])
toplabels.extend(['','Beads','','','','Lysis',''])

#### generate color table
colors=[color_one]
colors.extend(['white' for i in range(4)])
if n_beads == 3:
    colors.extend([color_two for i in range(n_beads)])
else:
    colors.extend([color_two for i in range(n_beads)])
    colors.extend(['white' for i in range(3-n_beads)])
colors.extend(['white'])

if n_lysis == 3:
    colors.extend([color_one for i in range(n_lysis)])
else:
    colors.extend([color_one for i in range(n_lysis)])
    colors.extend(['white' for i in range(3-n_lysis)])
########

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks([])
axT=ax.twiny()

barplot=ax.bar(x,y)
plt.ylim([-0.5,10.5])

for i in range(len(barplot)):
    barplot[i].set_color(colors[i])
    if colors[i]=='white':
        barplot[i].set_ec('k')

axT.set_xlim(ax.get_xlim())
ax.set_xticks(x)
axT.set_xticks(x)
axT.set_xticklabels(toplabels)

plt.show()
