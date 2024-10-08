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
press_list=['100', '1000', '10000']
structure=struct_list[3]
temperature=temp_list[1]
#pressure=press_list[1]

print(gas1,gas2, gas3)

fig = plt.figure(figsize=(10.8, 4.5))
fig.subplots_adjust(wspace=0.2)

for i, press in enumerate(press_list):

    file_name="data/data_IAST_%s-%s-%s_selectivity_%s_at_%s_P_%s.dat" %(gas1,gas2,gas3,structure,temperature,press)
    ax = fig.add_subplot(1, 3, i + 1, projection='ternary',ternary_sum=1)
#ax = plt.subplot(projection="ternary", ternary_sum=1)


    points=pd.read_csv(file_name, sep='\t', skiprows=1, header=None)
#    prenom=points[5]/points[6]
#    denom=points[2]/points[3]
    prenom=points[4]/points[5]
    denom=points[1]/points[2]
    selectivity=prenom/denom
    vmin = selectivity.min()
    vmax = selectivity.max()
#    vmin = 60
#    vmax = 160
    levels = np.linspace(vmin, vmax, 12)

    mode = 'axis'
    position = 'tick1'

    cmap = "PRGn"
    colors='r'
    shading = "gouraud"
    cs = ax.tricontour(points[1],points[2],points[3], selectivity,  cmap=cmap,levels=levels, linewidths=2.)

    #cs = ax.tripcolor(points[1],points[2],points[3], selectivity, cmap=cmap, shading=shading, rasterized=True)
    #ax.tricontour(points[1],points[2],points[3], points[4], color=color_key[gas3], linewidths=0.5)
    ax.clabel(cs)


    ax.set_tlabel(gas1)
    ax.set_llabel(gas2)
    ax.set_rlabel(gas3)
    ax.laxis.set_label_rotation_mode(mode)
    ax.raxis.set_label_rotation_mode(mode)
    ax.taxis.set_label_position(position)
    ax.laxis.set_label_position(position)
    ax.raxis.set_label_position(position)
    ax.legend()
    ax.taxis.set_minor_locator(MultipleLocator(10))
    ax.raxis.set_minor_locator(MultipleLocator(10))
    ax.laxis.set_minor_locator(MultipleLocator(10))
#ax.grid()
    ax.grid(axis='l', which='both', linestyle='--')
    ax.grid(axis='t', which='both', linestyle='--')
    ax.grid(axis='r', which='both', linestyle='--')
    alpha=05.3
#ax.set_title("${\\mathbf{\\alpha}}$ = " + str(alpha))
    pr=int(press)/1000
    #ax.set_title(" %s / %s \n P = %s bar" %(gas2, gas3,pr))
    ax.set_title(" %s / %s \n P = %s bar" %(gas1, gas2,pr))
#ax.grid(multiple=5, color="blue")


#cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9], transform=ax.transAxes)
#colorbar = plt.colorbar(cs, cax=cax)
#colorbar.set_label("PDF", rotation=270, va="baseline")

plt.savefig('Gas_%s_%s_selectivity_%s_at_%s.png' %(gas1,gas2,structure,temperature), bbox_inches='tight')
#plt.savefig('Gas_%s_%s_selectivity_%s_at_%s.png' %(gas2,gas3,structure,temperature), bbox_inches='tight')
plt.show()
