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
import sys
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
#press=press_list[0]

print(gas1,gas2, gas3)

x0=283
xn=313
nT=16
Textra=np.linspace(x0,xn,nT)
num_plots=nT


icol=0  #a random temperature step
temperature=Textra[icol]
temperature=temp_list[icol]
#selCO2_N2=pd.DataFrame()

#selCO2_N2=pd.DataFrame(columns=[['T', 'O2_rich', 'N2_rich', 'CO2_rich', 'equim']])
#selO2_N2=pd.DataFrame(columns=[['T','O2_rich','N2_rich','CO2_rich','equim']])

for i, press in enumerate(press_list):
#     for icol in range(len(Textra)):

#        fig = plt.figure(figsize=(10.8, 4.5))
#        fig.subplots_adjust(wspace=0.2)
        iloc1=0
        iloc2=42
        iloc3=945
        iloc4=525



#    for i, press in enumerate(press_list):
        file_name="data/data_IAST_%s-%s-%s_selectivity_%s_at_%s_P_%s.dat" %(gas1,gas2,gas3,structure,temp_list[icol],press)
#        ax = fig.add_subplot(1, 3, i + 1, projection='ternary',ternary_sum=1)
#ax = plt.subplot(projection="ternary", ternary_sum=1)


        points=pd.read_csv(file_name, sep='\t', skiprows=1, header=None)
        
        
        bulkP1a=points.loc[iloc1,1]/points.loc[iloc1,2]
        bulkP1b=points.loc[iloc1,3]/points.loc[iloc1,2]
        bulkP2a=points.loc[iloc2,1]/points.loc[iloc2,2]
        bulkP2b=points.loc[iloc2,3]/points.loc[iloc2,2]
        bulkP3a=points.loc[iloc3,1]/points.loc[iloc3,2]
        bulkP3b=points.loc[iloc3,3]/points.loc[iloc3,2]
        bulkP4a=points.loc[iloc4,1]/points.loc[iloc4,2]
        bulkP4b=points.loc[iloc4,3]/points.loc[iloc4,2]
        
        
        data=[{'T':Textra[icol],'O2_rich':points.loc[iloc1,4]/points.loc[iloc1,5]/bulkP1a,'N2_rich':points.loc[iloc2,4]/points.loc[iloc2,5]/bulkP2a,'CO2_rich':points.loc[iloc3,4]/points.loc[iloc3,5]/bulkP3a,'equim':points.loc[iloc4,4]/points.loc[iloc4,5]/bulkP4a}]
        
        
        data2=[{'T':Textra[icol],'O2_rich':points.loc[iloc1,6]/points.loc[iloc1,5]/bulkP1b,'N2_rich':points.loc[iloc2,6]/points.loc[iloc2,5]/bulkP2b,'CO2_rich':points.loc[iloc3,6]/points.loc[iloc3,5]/bulkP3b,'equim':points.loc[iloc4,6]/points.loc[iloc4,5]/bulkP4b}]
#        selCO2_N2=pd.DataFrame(data,index=[icol])
        if (i) == 0: 
            print(data)
            selCO2N2=pd.DataFrame(data)
            selO2N2=pd.DataFrame(data2)
#            print(selCO2N2)
#           selO2_N2=pd.DataFrame(data2)
        if (i) > 0:
            print(data)
            temp=pd.DataFrame(data)
            temp2=pd.DataFrame(data2)
            selCO2N2=pd.concat([selCO2N2,temp],ignore_index=True)
            selO2N2=pd.concat([selO2N2,temp2],ignore_index=True)
#            print(selCO2N2)
#           selCO2_N2=selCO2_N2.append(data, ignore_index=True)
#           selCO2_N2.loc[len(selCO2_N2)]=data
#           selO2_N2.loc[icol]=data2
#        selCO2_N2[['O2_rich']]=points.loc[0,4]/points.loc[0,5]/bulkP1a
#        selO2_N2[['O2_rich']]=points.loc[0,6]/points.loc[0,5]/bulkP1b
print("**********************************")
print(selCO2N2)
print("**********************************")
print(selO2N2)
selCO2N2_transposed=selCO2N2.T
selO2N2_transposed=selO2N2.T
file_out1="points1234/IAST_%s-%s_selectivity_%s_at_P_%s.dat" %(gas1,gas2,structure,press)
selCO2N2_transposed.to_csv(file_out1, sep='\t')
file_out2="points1234/IAST_%s-%s_selectivity_%s_at_P_%s.dat" %(gas3,gas2,structure,press)
selO2N2_transposed.to_csv(file_out2, sep='\t')
#print(data)
#print(selO2_N2)        
        #print(selCO2_N2['T'],selCO2_N2['O2_rich'],selO2_N2['O2_rich'])

sys.exit()

        

