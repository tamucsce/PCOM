
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
from pylab import *


font = {'family' : 'normal',
        'size'   : 40}



N = 5
prover_time_naive = (	0.048,	0.140625,	0.5625,	2.42188, 9.79688, 39.7969,	165.922,	667.703,	2661.91,	10688.7 , 42754.8, 171019.2,684076.8)
verfier_time_naive = (	0.02,	0.03125,	0.09375,	0.204002,	0.375,	0.8125,	1.82812,	3.92188,	8.5,	18.2031 ,38.1094,76.21,152.42)
prover_time_fft = (	0.421875,	1.67188,	5.5625,	17.7969,	52.625,	150.891,	399.391,	1068.47,	2704.36, 6844.2 , 17740.8, 43307.8, 108269)
verfier_time_fft = (	0.015625,	 0.04,	0.078125,	0.171875,	0.390625,	0.84375,	1.6875,	3.79688,	8.1875,	17.375 ,38.1094,76.21,152.42)



proof_size = (12288,	30720,	73728,	172032,	393216,	884736,	1966080,	4325376,	9437184,	20447232 ,44040192,95693045.7, 206696978)
proof_size = (57344,
143360,
344064,
802816,
1835008,
4128768,
9175040,
20185088,
44040192,
95420416,
205520896,
440401920,
939524096,
)

x = (2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16)
#x1=(10)
#x2=(100,1000,10000,100000,200000)

fig, ax = plt.subplots(figsize=(23,15))
ax2 = ax.twinx()

rects1 = ax.plot(x, prover_time_naive,color='k',linewidth=3,marker='o',markersize=20,fillstyle='full')
rects2 = ax.plot(x, verfier_time_naive,color='b',linewidth=3,marker='^',markersize=20,fillstyle='full')
rects4 = ax.plot(x, prover_time_fft,color='g',linewidth=3,marker='d',markersize=20,fillstyle='full')


ax.set_yscale('log')
ax.set_xscale('log')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=30)
ax.xaxis.set_ticks((2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384))
ax.set_xticklabels( ('$2^1$', '$2^2$', '$2^3$', '$2^4$', '$2^5$','$2^6$', '$2^7$', '$2^8$', '$2^9$', '$2^{10}$', '$2^{11}$', '$2^{12}$', '$2^{13}$') )

ax.xaxis.set_ticks((0,0,16,64,256,1024,4096,16384, 65536))
ax.set_xticklabels( ('' , '', '$2^4$', '$2^6$',  '$2^8$', '$2^{10}$',  '$2^{12}$', '$2^{14}$', '$2^{16}$',  ) )

ax.yaxis.set_ticks((0.01,0.01,0.01,1,100,10000,1000000,))
ax.set_yticklabels(('' , '', '$10^{-2}$', '$10^0$', '$10^2$', '$10^4$', '$10^6$','$10^8$', '$10^{10}$'),ha='left')
plt.rcParams.update({'legend.labelspacing':0.25})


ax2.set_yscale('log')
ax2.yaxis.set_ticks((0.01,0.01,10000,100000,1000000,10000000,100000000,1000000000,10000000000, 100000000000, 1000000000000,1000000000000,10000000000000,100000000000000,1000000000000000))
#ax2.set_ylim([0,10000])

ax2.set_ylim([45000,1000000000000000])
ax2.set_yticklabels(('','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$','$10^6$','$10^7$','$10^8$','$10^9$','$10^{10}$','$10^{10}$','$10^{11}$','$10^{12}$','$10^{13}$'),ha='left',fontsize=50)

rects3 = ax2.plot(x, proof_size,color='r',linewidth=3,marker='s',markersize=20,fillstyle='full')
leg = ax.legend( (rects1[0], rects4[0], rects2[0], rects3[0]), ('Prover Time Naive','Prover Time MultiPoint', 'Verifier Time','Proof Size') ,loc='upper left', borderpad=0.1,bbox_to_anchor=[0, 1], fontsize=55)

ax.yaxis.grid(True)
ax.xaxis.grid(True)




ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='y', pad=90)

#ax.text(150000, 6, '$2\\times10^5$', fontsize=50)


for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(50) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
ax2.set_ylabel('Proof Size (kbytes)',fontsize=40,**font,labelpad=30)
ax.set_ylabel('Time (seconds)',fontsize=40,**font,labelpad=30)
ax.set_xlabel('degree of poly ',fontsize=40,**font,labelpad=10)


plt.show()
plt.savefig('plot.pdf',dpi=120)
