import numpy as np
import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

import pylab
from pylab import *






N = 4


ind = np.arange(N)  # the x locations for the groups
print(ind)
width = 0.11      # the width of the bars

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.125, bottom=0.16, right=0.9, top=0.995,wspace=0, hspace=0)




prover = [293631,294458, 294263, 293388, ]
rects2 = ax.bar(ind, prover, width, color='blue', log =False,linewidth=0.5,)

multipoint_eval = [302893, 303346, 302739, 301973, ]
rects3 = ax.bar(ind, multipoint_eval, width, color='limegreen', log =False,linewidth=0.5, bottom = prover)

lamda = [643,644,1144, 1316,]
rects1 = ax.bar(ind, lamda, width, color='black',log =False,linewidth=0.5,  bottom = multipoint_eval)

decryption = [2484, 2637, 3002, 4151,]
rects4 = ax.bar(ind, decryption, width, color='red', log =False,linewidth=0.5, bottom = lamda)

#prover_all = np.add(lamda, prover,multipoint_eval,decryption).tolist()

# verfier time and bars
verifier = [2812,2799,2959,3050, ]
rects5 = ax.bar(ind, verifier, width, color='orange', log =False,linewidth=0.1, bottom =decryption)

Encryption = [2484 ,2637,3002,4151]
rects6 = ax.bar(ind, Encryption, width, color='purple', log =False,linewidth=0.1, bottom = verifier)

# planar1 = (1.353008174,7.50964737,13.82530701,59.07403015,513.2520905)
# rects5 = ax.bar(ind+2*width, planar1, width, color='yellow',edgecolor='yellow', log =True,linewidth=0.1)

# planar2 = (0.853008174,5.50964737,10.82530701,39.07403015,413.2520905)
# rects6 = ax.bar(ind+3*width, planar2, width, color='w', log =True,linewidth=2)

#ax.yaxis.set_ticks((0.01,1,100,10000,1000000,100000000,10000000000))
#ax.set_yticklabels(('$10^{-2}$', '$10^0$', '$10^2$', '$10^4$', '$10^6$','$10^8$', '$10^{10}$'),ha='left')

ax.yaxis.set_ticks((0,0,100000,200000,300000,400000,500000,600000,))
ax.set_yticklabels(('','0','$100$', '$200$', '$300$', '$400$', '$500$','$600$','$700$','$800$', ),ha='left',)
# plt.rcParams.update({'legend.labelspacing':0.25})

# # add some
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=30)
#ax.set_title('Scores by group and gender')
ax.set_xticks(ind+width)
ax.xaxis.set_ticks((0.0,0.0,0,1,2,3,))
ax.set_xticklabels( ('','','$2$', '$8$', '$32$', '$128$', ) )

#ax.set_yscale('log')
#leg = ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('Semihonest', 'Polynomial Commitment','\\textsc{Buffet} BFS','Certifying Algorithm','Planar Separator 1', 'Planar Separator 2') ,bbox_to_anchor=(0., 1.005, 1., .1005),loc=3,ncol=2,mode="expand",borderaxespad=0,fontsize = 30)
leg = ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('Lambda', 'Prover', 'Multipoint Eval' ,'Decryption', 'verifier','Encryption') ,bbox_to_anchor=(0., 1.005, 1., .1005),loc=3,ncol=2,mode="expand",borderaxespad=0,fontsize = 40)
# leg = ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('Strawman', '\\textsc{Libsnark} BFS','\\textsc{Buffet} BFS','Certifying Algorithm','Planar Separator 1', 'Planar Separator 2') ,bbox_to_anchor=(0., 1.005, 1., .1005),loc=3,ncol=2,mode="expand",borderaxespad=0,fontsize = 40)

ax.yaxis.grid(True)

# rects2[2].set_color('w')
# rects2[2].set_edgecolor('b')
# rects2[2].set_hatch('//')
# rects2[2].set_linewidth(2)

# rects2[3].set_color('w')
# rects2[3].set_edgecolor('b')
# rects2[3].set_hatch('//')
# rects2[3].set_linewidth(2)

# #rects2[4].set_color('w')
# #rects2[4].set_edgecolor('b')
# #rects2[4].set_hatch('//')
# #rects2[4].set_linewidth(2)

# rects2[1].set_color('w')
# rects2[1].set_edgecolor('b')
# rects2[1].set_hatch('//')
# rects2[1].set_linewidth(2)

# #rects1[4].set_color('w')
# #rects1[4].set_edgecolor('k')
# #rects1[4].set_hatch('//')
# #rects1[4].set_linewidth(2)

# #rects4[4].set_color('w')
# #rects4[4].set_edgecolor('#00FFFF')
# #rects4[4].set_hatch('//')
# #rects4[4].set_linewidth(2)


# rects3[2].set_color('w')
# rects3[2].set_edgecolor('limegreen')
# rects3[2].set_hatch('//')
# rects3[2].set_linewidth(2)

# rects3[3].set_color('w')
# rects3[3].set_edgecolor('limegreen')
# rects3[3].set_hatch('//')
# rects3[3].set_linewidth(2)

#rects3[4].set_color('w')
#rects3[4].set_edgecolor('limegreen')
#rects3[4].set_hatch('//')
#rects3[4].set_linewidth(2)
#patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')
#for bar, pattern in zip(rects2, patterns):
#     bar.set_hatch(pattern)
#def autolabel(rects):
#    # attach some text labels
#    for rect in rects:
#        height = rect.get_height()
#        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
#                ha='center', va='bottom')
#
#autolabel(rects1)
#autolabel(rects2)
plt.xlim((-0.2, 4.5))
plt.ylim((100,650000))
ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='y', pad=100)




for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(50) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
ax.set_ylabel('Time (s)',fontsize=50)
ax.set_xlabel('Num Parties',fontsize=50)

plt.subplots_adjust(left=0.13, bottom=0.17, top=0.77, right=0.98)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
ax.plot()
fig.set_size_inches(20,10)
plt.savefig('bar.png')

