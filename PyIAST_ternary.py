# coding: utf-8

# # Test

# We test the IAST code with a binary mixture of methane (CH$_4$) and ethane (C$_2$H$_6$) adsorbed in metal-organic framework IRMOF-1.
#
# Simulated pure-component adsorption isotherms at 298 K are present in the files:
# * `IRMOF-1_methane_isotherm_298K.csv`
# * `IRMOF-1_ethane_isotherm_298K.csv`
#
# We ran dual component GCMC mixture isotherms of methane/ethane in IRMOF-1 at 65 bar total pressure and 298 K at different mixture compositions. This data is present in `IRMOF-1_methane_ethane_mixture_isotherm_65bar_298K.csv`.
#
# The goal of this test is to use pyIAST to predict the mixture isotherms from the pure-component isotherms and compare to the binary GCMC mixture simulations.

# In[1]:

from __future__ import absolute_import
from __future__ import print_function
import pyiast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mpltern
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import sys
from six.moves import range

import matplotlib.tri as tri
# Matplotlib settings
plt.style.use('bmh')

# colors to use for plots
color_key = {'CO2': 'g', 'N2': 'b', 'O2': 'r'}
struct_list = ['Alumina', 'WS2050', 'XLBFK-12', 'XLBFK-16']
temp_list=['283K', '298K', '313K']
gas_list=['CO2', 'N2', 'O2']
gas1=gas_list[0]
gas2=gas_list[1]
gas3=gas_list[2]

structure=struct_list[3]
temperature=temp_list[0]

print(gas1,gas2, gas3)
# ## Load pure-component isotherm data as Pandas dataframes

# In[2]:

df_gas1 = pd.read_csv("isotherms/%s_%s_isotherm_%s.dat" %(structure, gas1, temperature), sep='\s+')
df_gas1.head()

# In[3]:

df_gas2 = pd.read_csv("isotherms/%s_%s_isotherm_%s.dat" %(structure, gas2, temperature), sep='\s+')
df_gas2.head()

df_gas3 = pd.read_csv("isotherms/%s_%s_isotherm_%s.dat" %(structure, gas3, temperature), sep='\s+')
df_gas3.head()


num=80
dp=0.03
#total pressure in bar
#total_pressure=1000
mol_frac=np.linspace(dp,1-dp,num)
# ### Plot isotherm data

a=pd.DataFrame(mol_frac, columns=['%s_mole_fract' %gas2])
b=pd.DataFrame(mol_frac, columns=['%s_mole_fract' %gas2])
         
selectivity=pd.DataFrame(mol_frac, columns=['%s_mole_fract' %gas2])
bulk =(1-a['%s_mole_fract' %gas2])/a['%s_mole_fract' %gas2]
# In[4]:

fig, ax = plt.subplots()

plt.scatter(
    df_gas2['p(mbar)'],
    df_gas2['loading(mmol/g)'],
    label=gas2,
    color=color_key[gas2],
    s=50)
plt.scatter(
    df_gas1['p(mbar)'],
    df_gas1['loading(mmol/g)'],
    label=gas1,
    color=color_key[gas1],
    s=50,
    marker='s')
    
plt.scatter(
    df_gas3['p(mbar)'],
    df_gas3['loading(mmol/g)'],
    label=gas3,
    color=color_key[gas3],
    s=50,
    marker='s')

plt.xlabel('Pressure (mbar)')
plt.ylabel('Gas uptake (loading(mmol/g))')

plt.xlim([0, 20000])
plt.ylim(ymin=0)

plt.legend(loc='lower right')

plt.tight_layout()
#plt.savefig("pure_component_isotherms.pdf", format='pdf')
plt.savefig("graphs/pure_component_%s_and_%s_isotherms_of_%s_at_%s.png"  
           %(gas1, gas2, structure, temperature), format='png', dpi=250)
plt.show()
###############

# ### Interpolator isotherm for Methane

# In[7]:
    
gas1_isotherm = pyiast.ModelIsotherm(
#gas1_isotherm = pyiast.InterpolatorIsotherm(
    df_gas1,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
#    fill_value=df_gas1["loading(mmol/g)"].max())
#    model="Langmuir")
#    model="TemkinApprox", param_guess= {'M':10, 'K':0.03,'theta':0.0001})
#(0)         
#    model="BET", param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
#(1) (default)          
#    model="BET")
#(2)                          
#    model="DSLangmuir", param_guess= {'M1':3.36, 'K1':0.002,'M2':3.2,'K2':1.78})
#(3)                  
#     model="DSLangmuir", param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#(3)                 
    model="BET", param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
    #, param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#    fill_value=df_gas1['p(mbar)'].max())

param_names = [param for param in gas1_isotherm.params.keys()]
print ('param_names for %s:' %gas1, param_names)
# ### Interpolator isotherm for ethane
gas1_isotherm.print_params()
pyiast.plot_isotherm(gas1_isotherm)
# In[8]:
####******************************************************************************************************66
#sys.exit()
####******************************************************************************************************77
gas2_isotherm = pyiast.ModelIsotherm(
#gas2_isotherm = pyiast.InterpolatorIsotherm(
    df_gas2,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
#    fill_value=df_gas2["loading(mmol/g)"].max())
#    model="Langmuir")
    model="Langmuir")
    #model="DSLangmuir", param_guess= {'M1':1300, 'K1':10,'M2':1000,'K2':10})
#    fill_value=df_gas2["p(mbar)"].max())


param_names = [param for param in gas2_isotherm.params.keys()]
print ('param_names for %s:' %gas2, param_names)

gas2_isotherm.print_params()

pyiast.plot_isotherm(gas2_isotherm)

gas3_isotherm = pyiast.ModelIsotherm(
#gas2_isotherm = pyiast.InterpolatorIsotherm(
    df_gas3,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
    model="Langmuir")
#    model="DSLangmuir")


param_names = [param for param in gas3_isotherm.params.keys()]
print ('param_names for %s:' %gas3, param_names)

gas3_isotherm.print_params()
pyiast.plot_isotherm(gas3_isotherm)



#partial_pressures=[0.166, 0.679, 0.02]
#component_loadings=pyiast.iast(partial_pressures, [gas1_isotherm, gas2_isotherm, gas3_isotherm])

#print ('print component_loadings:' , component_loadings)

ints=45
intsq=(ints-2)*(ints-3)
mol_fract = np.zeros((3, intsq))
index=0
for id_ in range(1,ints):
  
  for ii in range(1,ints):
      if   (ints-id_-ii) > 0: 
           mol_fract[0,index]=id_
           mol_fract[1,index]=ii
           mol_fract[2,index]=ints-id_-ii
           index=index+1
#mol_fract=mol_fract/(ints)
print ("Original dataset")
print (mol_fract)      

mol_fract=mol_fract[:,~np.all(mol_fract == 0, axis=0)]
#mol_fract=mol_fract[~np.all(mol_fract == 0, axis=1)]
print ("After removing rows")

mol_fract=mol_fract/(ints)
print (mol_fract)
#sys.exit()
#sizea=mol_fract.size/3
#mol_fract = np.ma.masked_equal(mol_fract,0)
#mol_fract=mol_fract.compressed()
#print (mol_fract.size)

#sizeb=mol_fract.size/3
#print ('number of rows with zero',sizea-sizeb)
#sys.exit()
#########
size=int(mol_fract.size/3)
print (size)

###########                     pressure is in mbar #########################################
press_list=[100,1000,10000]
total_pressure=press_list[2]
partial_pressures=mol_fract*total_pressure
print(partial_pressures)
#sys.exit()
#mbar=1000*total_pressure
#a=pd.DataFrame(mol_frac, columns=['%s_mole_fract' %gas2])
iast_component_loadings = np.zeros((3,size))  # store component loadings here
#iast_component_loadings = pd.DataFrame(columns=['%s_mole_fract' %gas1,'%s_mole_fract' %gas2,'%s_mole_fract' %gas3])
# np.zeros((3, intsq-1))  # store component loadings here

#component_loadings=pyiast.iast(partial_pressures, [gas1_isotherm, gas2_isotherm, gas3_isotherm])

for i in range(size):
              #y_gas2 = mol_frac[i]
#              partial_pressures = total_pressure*np.array([y_gas2, 1.0 - y_gas2])
              iast_component_loadings[:, i] = pyiast.iast(
              partial_pressures[:,i], [gas1_isotherm, gas2_isotherm, gas3_isotherm], verboseflag=False)
              print ('print component_loadings at index %s:' %i, iast_component_loadings[:,i])

file_name="data/data_IAST_%s-%s-%s_selectivity_%s_at_%s_P_%s.dat" %(gas1,gas2,gas3,structure,temperature,total_pressure)
#a=pd.DataFrame(Loading1, columns=['#', 'ZIF69'])
#DF = pd.DataFrame(iast_component_loadings)
columns = ['x_%s' %gas1,'x_%s' %gas2,'x_%s' %gas3]
DF=pd.DataFrame(mol_fract.T, columns=columns)
newcolumns = ['ads_%s' %gas1,'ads_%s' %gas2,'ads_%s' %gas3]
DF2=pd.DataFrame(iast_component_loadings.T, columns=newcolumns)
DF=DF.assign(**DF2)
DF.to_csv(file_name, sep='\t')


#T=tri.Triangulation(DF['x_%s' %gas1],DF['x_%s' %gas2])
#plt.tricontourf(DF['x_%s' %gas1],DF['x_%s' %gas2],DF['x_%s' %gas3],DF['ads_%s' %gas1])

#plt.show()
#sys.exit()
#iast_component_loadings[~np.isnan(iast_component_loadings).any(axis=1)]
#print ('print component_loadings', a[:,:])






# In[16]:

pyiast._MODELS
