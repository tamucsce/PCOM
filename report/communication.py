
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
from pylab import *


font = {'family' : 'normal',
        'size'   : 40}



N = 5
com_2party = (24.048, 42.568, 77.520, 145.288, 278.640, 543.112, 6979.120, 12713.280, 25891.680, 52481.040,105892.320,211707.120, 425352.240, 861944.880, 1746370.560, 3538167.840,7076335.680)
com_128party = (3078.144, 5449, 9922.560,18596.864 ,35665.920, 69518,803727,  1627299, 3314135, 6717573,13554216, 27098511, 54445086,109400172,218800344,437600688, 875201376)
com_500party = (12024, 38760, 72644, 139320, 271556, 3139560,6356640, 12945840,26240520, 52946160, 105853560,211707120,425352240,850704480, 1701408960,3402817920,6805635840)

com_2party = (0.024048,  0.042568, 0.07752, 0.145288, 0.27864, 0.543112,        6.27912, 12.71328  ,      25.89168,        52.48104,        105.89232   ,    211.70712 ,      425.35224  )
com_128party = (3.054096  , 5.406136,  9.84504, 18.451576 , 35.38728 ,68.975224,  803.72736,  1627.2998, 3314.13504, 6717.573,  13554.2169,  27098.51, 54445.086,  )
com_500party = (11.7187031, 20.743585, 37.7758593,   70.799523,  135.7825, 264.66102,  3065.9765, 6207.656,   12642.421,  25625.5078,  51705.234,  103372.61, 207691.523)
com_1000party = (23.46,   41.528,   75.62,    141.74,    271.84,  529.85,   6131.95,    12415.3,    25284.84,    51251.01,   103410.47,   206745.23,     415383.05)


x=(10,100,1000,10000,100000,200000)
x = (2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16,2**17,2**18,2**19,2**20)
x = (2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16)

#x1=(10)
#x2=(100,1000,10000,100000,200000)

fig, ax = plt.subplots(figsize=(22,15))
#ax2 = ax.twinx()

rects1 = ax.plot(x, com_2party,color='k',linewidth=3,marker='o',markersize=20,fillstyle='full')
rects2 = ax.plot(x, com_128party,color='b',linewidth=3,marker='^',markersize=20,fillstyle='full')
rects3 = ax.plot(x, com_1000party,color='g',linewidth=3,marker='d',markersize=20,fillstyle='full')

#plt.scatter(x1,Dynamic1, s=500, marker='o',facecolor='k')
#plt.scatter(x2,Dynamic2, s=500, marker='o',edgecolor='k',linewidth='3', facecolor='w', hatch='////')


ax.set_yscale('log')
ax.set_xscale('log')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=30)
ax.xaxis.set_ticks((2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576))
ax.set_xticklabels( ('$2^1$', '$2^2$', '$2^3$', '$2^4$', '$2^5$','$2^6$', '$2^7$', '$2^8$', '$2^9$', '$2^{10}$', '$2^{11}$', '$2^{12}$', '$2^{13}$', '$2^{14}$', '$2^{15}$','$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$',) )

ax.xaxis.set_ticks((0,0,4,16,64,256,1024,4096,16384,65536))
ax.set_xticklabels( ('' , '', '$2^2$', '$2^4$',  '$2^6$', '$2^8$',  '$2^{10}$',  '$2^{12}$', '$2^{14}$', '$2^{16}$') )

ax.yaxis.set_ticks((0.001 , 0.001,0.01,1,100,10000,1000000,100000000,10000000000,))
ax.set_yticklabels(('' , '','$10^{-2}$','$10^0$', '$10^2$','$10^4$','$10^6$', '$10^8$', '$10^{10}$'),ha='left')
plt.rcParams.update({'legend.labelspacing':0.25})


#ax2 = ax.twinx()
# ax2.yaxis.set_ticks((0,0,1000,5000,10000,15000,20000,25000,30000,35000,40000,45000))
# ax2.set_ylim([0,45000])
# ax2.set_yticklabels(('','','$1$', '$5$','${10}$','${15}$','${20}$','${25}$','${30}$','${35}$','${40}$','${45}$',),ha='left',fontsize=50)

leg = ax.legend( (rects1[0],rects2[0],rects3[0]), ('2 parties','128 parties','1000 parties',) ,loc='upper left', borderpad=0.1,bbox_to_anchor=[0, 1],fontsize=55)

ax.yaxis.grid(True)
ax.xaxis.grid(True)


plt.ylim((0.01, 1000000))

plt.xlim((10, 80000))

ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='y', pad=90)

#ax.text(150000, 6, '$2\\times10^5$', fontsize=50)


for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(50) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(50)
# for tick in ax2.yaxis.get_major_ticks():
#                 tick.label.set_fontsize(50)
#ax2.set_ylabel('Proof Size (kbytes)',fontsize=60,**font,labelpad=30)
ax.set_ylabel('Communication cost (MB)',fontsize=60,**font,labelpad=30)
ax.set_xlabel('degree of poly',fontsize=60,**font,labelpad=10)


plt.show()
plt.savefig('plot.pdf',dpi=120)
