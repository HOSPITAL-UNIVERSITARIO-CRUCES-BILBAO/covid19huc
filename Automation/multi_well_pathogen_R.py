import matplotlib.pyplot as plt

###################
wash_two=4
wash_one=3
color_one='yellowgreen'
color_two='indianred'
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
if wash_two == 3:
    colors.extend([color_two for i in range(wash_two)])
else:
    colors.extend([color_two for i in range(wash_two)])
    colors.extend(['white' for i in range(4-wash_two)])
colors.extend(['white'])

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
