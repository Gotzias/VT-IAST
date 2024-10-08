"""
======
Legend
======

The legend can be put using ``ax.legend`` in the same way as Matplotlib.
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import mpltern
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# colors to use for plots
color_key = {'CO2': 'g', 'N2': 'b', 'O2': 'r'}
struct_list = ['Alumina', 'WS2050', 'XLBFK-12', 'XLBFK-16']
temp_list=['283K', '298K', '313K']
gas_list=['CO2', 'N2', 'O2']
gas1=gas_list[0]
gas2=gas_list[1]
gas3=gas_list[2]
#gasa is either CO2 or O2
gasa=gas_list[2]
#gasb is N2 (always)
gasb=gas_list[1]
press_list=['100', '1000', '10000']
structure=struct_list[3]
temperature=temp_list[0]
#pressure=press_list[1]

print(gas1,gas2, gas3)

fig = plt.figure(figsize=(10, 4))
fig.subplots_adjust(wspace=0.2)

for i, press in enumerate(press_list):

    file_name="data/data_IAST_%s-%s-%s_selectivity_%s_at_%s_P_%s.dat" %(gas1,gas2,gas3,structure,temperature,press)
    ax = fig.add_subplot(1, 3, i + 1, projection='ternary',ternary_sum=1)
#ax = plt.subplot(projection="ternary", ternary_sum=1)


    points=pd.read_csv(file_name, sep='\t', skiprows=1, header=None)

    prenom=points[4]/points[5]
    denom=points[1]/points[2]
    if gasa==gas3:
        prenom=points[6]/points[5]
        denom=points[3]/points[2]
    selectivity=prenom/denom
    vmin = selectivity.min()
    vmax = selectivity.max()
    print(vmin,np.log10(vmin))
    print(vmax,np.log10(vmax))
    vmin = 0.26
    vmax = 1.56
#    levels = np.linspace(vmin, vmax, 7)
    levels = np.logspace(np.log10(vmin), np.log10(vmax), num=10, endpoint=True,base=10.0) 

    mode = 'axis'
    position = 'tick1'

    cmap = "PRGn"
    colors='black'
    shading = "gouraud"

#    cs = ax.tricontour(points[1],points[2],points[3], selectivity,  cmap=cmap,levels=levels, linewidths=2.)
    cs = ax.tricontourf(points[1],points[2],points[3], selectivity,  cmap=cmap,levels=levels)
#    cs = ax.tripcolor(points[1],points[2],points[3], selectivity,  shading=shading, cmap=cmap)

    #cs = ax.tripcolor(points[1],points[2],points[3], selectivity, cmap=cmap, shading=shading, rasterized=True)
    #ax.tricontour(points[1],points[2],points[3], points[4], color=color_key[gas3], linewidths=0.5)
    ax.clabel(cs,colors=colors,fontsize='smaller',inline=False)


    ax.set_tlabel(gas1)
    ax.set_llabel(gas2)
    ax.set_rlabel(gas3)
    ax.laxis.set_label_rotation_mode(mode)
    ax.raxis.set_label_rotation_mode(mode)
    ax.taxis.set_label_position(position)
    ax.laxis.set_label_position(position)
    ax.raxis.set_label_position(position)
#    ax.legend('Gotzias')
    ax.taxis.set_minor_locator(MultipleLocator(10))
    ax.raxis.set_minor_locator(MultipleLocator(10))
    ax.laxis.set_minor_locator(MultipleLocator(10))
#ax.grid()
    ax.grid(axis='l', which='both', linestyle='--')
    ax.grid(axis='t', which='both', linestyle='--')
    ax.grid(axis='r', which='both', linestyle='--')

    pr=int(press)/1000
    ax.set_title(" %s / %s selectivity\n P = %s bar" %(gasa, gasb,pr))
    if i==2:
       cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9], transform=ax.transAxes)
       colorbar = fig.colorbar(cs, cax=cax)
       colorbar.set_label(" %s / %s selectivity" %(gasa, gasb), rotation=270, va='baseline')
#ax.grid(multiple=5, color="blue")


#cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9], transform=ax.transAxes)
#colorbar = plt.colorbar(cs, cax=cax)
#colorbar.set_label("PDF", rotation=270, va="baseline")

plt.savefig('Gas_%s_%s_selectivity_%s_at_%s.png' %(gasa,gasb,structure,temperature), bbox_inches='tight')
plt.show()
