
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
from pylab import *


font = {'family' : 'normal',
        'size'   : 40}



N = 5
prover_time = (	0.1764,  0.3376 , 0.649 ,  1.301 ,  2.6  ,   5.26 ,   11.0404 ,22.040784   ,    44.00851264 ,    87.88475005 ,    175.5323946   ,  350.645035  ,    700.5586461)
#Dynamic1 = (2.112)
#Dynamic2 = (151.8744,16172.7744,1628711.774,162987101.8,651974201.8)
verfier_time = (0.017  ,   0.035  , 0.068   ,0.1312 , 0.247  , 0.484  , 0.9309 , 1.8273691 , 3.589, 7.054,   13.87,    27.298,   53.74)
proof_size = (768,	960,	1152,	1344,	1536,	1728,	1920,	2112,	2304,	2496,	2688,	2880,	3072,	3264,	3456,	3648,	3840)
proof_size = (4800,      5952,    7104,    8256,    9408,    10560,   11712,   12864,   14016,   15168,   16320,   17472,   18624)
x=(10,100,1000,10000,100000,200000)
x = (2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16,2**17,2**18,2**19,2**20)
x = (2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16)

#x1=(10)
#x2=(100,1000,10000,100000,200000)

fig, ax = plt.subplots(figsize=(22,15))
ax2 = ax.twinx()

rects1 = ax.plot(x, prover_time,color='k',linewidth=3,marker='o',markersize=20,fillstyle='full')
#plt.scatter(x1,Dynamic1, s=500, marker='o',facecolor='k')
#plt.scatter(x2,Dynamic2, s=500, marker='o',edgecolor='k',linewidth='3', facecolor='w', hatch='////')
rects2 = ax.plot(x, verfier_time,color='b',linewidth=3,marker='^',markersize=20,fillstyle='full')


ax.set_yscale('log')
ax.set_xscale('log')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=30)
ax.xaxis.set_ticks((2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576))
ax.set_xticklabels( ('$2^1$', '$2^2$', '$2^3$', '$2^4$', '$2^5$','$2^6$', '$2^7$', '$2^8$', '$2^9$', '$2^{10}$', '$2^{11}$', '$2^{12}$', '$2^{13}$', '$2^{14}$', '$2^{15}$','$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$',) )

ax.xaxis.set_ticks((0,0,4,16,64,256,1024,4096,16384,65536,262144,1048576))
ax.set_xticklabels( ('' , '', '$2^2$', '$2^4$',  '$2^6$', '$2^8$',  '$2^{10}$',  '$2^{12}$', '$2^{14}$', '$2^{16}$',) )

ax.yaxis.set_ticks((0.001,0.001, 0.01,0.1,1,10,100,1000,10000,))
ax.set_yticklabels(('' , '', '$10^{-2}$', '$10^{-1}$','$10^0$', '$10^1$','$10^2$','$10^3$', '$10^4$', '$10^6$','$10^8$', '$10^{10}$'),ha='left')
plt.rcParams.update({'legend.labelspacing':0.25})


#ax2 = ax.twinx()
ax2.yaxis.set_ticks((0,0,1000,5000,10000,15000,20000,25000,30000,35000,40000,45000))
ax2.set_ylim([0,45000])
ax2.set_yticklabels(('','','$1$', '$5$','${10}$','${15}$','${20}$','${25}$','${30}$','${35}$','${40}$','${45}$',),ha='left',fontsize=50)

rects3 = ax2.plot(x, proof_size,color='r',linewidth=3,marker='s',markersize=20)
leg = ax.legend( (rects1[0], rects2[0], rects3[0]), ('Prover Time','Verifier Time','Proof Size') ,loc='upper left', borderpad=0.1,bbox_to_anchor=[0, 1],fontsize=55)

ax.yaxis.grid(True)
ax.xaxis.grid(True)



plt.xlim((10, 80000))

ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='y', pad=90)

#ax.text(150000, 6, '$2\\times10^5$', fontsize=50)


for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(50) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
ax2.set_ylabel('Proof Size (KB)',fontsize=40,**font,labelpad=30)
ax.set_ylabel('Time (seconds)',fontsize=40,**font,labelpad=30)
ax.set_xlabel('degree of poly',fontsize=40,**font,labelpad=10)


plt.show()
plt.savefig('plot.pdf',dpi=120)
