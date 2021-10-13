# -*- coding: cp936 -*-
import matplotlib.pyplot as plt
import numpy as np
a = [293631,294458, 294263, 293388,   293611.6,  293519.2 ]
b = [302893, 303346, 302739, 301973,   301559.3, 301222.6 ]
c = [643,644,1144, 1316,     1818,2070]
d = [2484, 2637, 3002, 4151,   4946,5483 ]
e = [2812,2799,2959,3050,     2865,2986]
f = [2484 ,2637,3002,4151,    4946, 5483]
N=6
ind = np.arange(N)
b1 = list()
b2 = list()
b3 = list()
b4 = list()
for i in range (4):
    b1.append(a[i]+b[i])
    b2.append(a[i]+b[i]+e[i])
    b3.append(a[i]+b[i]+d[i]+e[i])
    b4.append(a[i]+b[i]+e[i]+d[i]+f[i])
    
#fig = plt.figure(figsize=(9,21))

fig, ax = plt.subplots()

plt.rcParams['savefig.dpi'] = 90 #图片像素
plt.rcParams['figure.dpi'] = 90
width = 0.36
plt.bar(ind, b,width,label='Evaluation')
plt.bar(ind, a,width,label='Proving',bottom=b)
plt.bar(ind, e,width,label='Verifying',bottom=b1,color='black')
plt.bar(ind, d,width,label='Decryption',bottom=b2,color='red')
plt.bar(ind, f,width,label='Encryption',bottom=b3,color='green')
plt.bar(ind, c,width,label='Random point',bottom=b4,color='purple')

plt.legend(bbox_to_anchor=(0., 1.005, 1., .1005),loc=3,ncol=3,mode="expand",borderaxespad=0,fontsize = 35)

ax.yaxis.set_ticks((0,0,100000,200000,300000,400000,500000,600000,))
ax.set_yticklabels(('','0','$100$', '$200$', '$300$', '$400$', '$500$','$600$','$700$','$800$', ),ha='left',)
# plt.rcParams.update({'legend.labelspacing':0.25})
#leg = ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('Lambda', 'Prover', 'Multipoint Eval' ,'Decryption', 'verifier','Encryption') ,bbox_to_anchor=(0., 1.005, 1., .1005),loc=3,ncol=2,mode="expand",borderaxespad=0,fontsize = 40)
# # add some
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif', size=30)
#ax.set_title('Scores by group and gender')
#ax.set_xticks(ind+width)
ax.xaxis.set_ticks((0.0,0.0,0,1,2,3,))
ax.set_xticklabels( ('','','$2$', '$8$', '$32$', '$128$', ) )

#ax.yaxis.grid(True)

plt.xlim((-0.5, 3.5))
plt.ylim((100,650000))
ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='y', pad=100)


for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(50) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
ax.set_ylabel('Time (s)',fontsize=50)
ax.set_xlabel('Numer of Parties',fontsize=50)

plt.subplots_adjust(left=0.17, bottom=0.17, top=0.85, right=0.95)
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
fig.set_size_inches(16,14)
plt.savefig('C:/Users/star/Desktop/mpsi/breakdown.pdf')