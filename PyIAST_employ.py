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
import sys
from six.moves import range

# Matplotlib settings
plt.style.use('bmh')

# colors to use for plots
color_key = {'CO2': 'g', 'N2': 'b', 'O2': 'r'}
struct_list = ['Alumina', 'WS2050', 'XLBFK-12', 'XLBFK-16']
temp_list=['283K', '298K', '313K']
gas_list=['CO2', 'N2', 'O2']
gas1=gas_list[0]
gas2=gas_list[1]

structure=struct_list[3]
temperature=temp_list[1]

print(gas1,gas2)
# ## Load pure-component isotherm data as Pandas dataframes

# In[2]:

df_gas1 = pd.read_csv("isotherms/%s_%s_isotherm_%s.dat" %(structure, gas1, temperature), sep='\s+')
df_gas1.head()

# ### Ethane

# In[3]:

df_gas2 = pd.read_csv("isotherms/%s_%s_isotherm_%s.dat" %(structure, gas2, temperature), sep='\s+')
df_gas2.head()


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
    model="BET", param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
#(1) (default)      
#    model="BET")
#(2)             
#    model="DSLangmuir", param_guess= {'M1':3.36, 'K1':0.002,'M2':3.2,'K2':1.78})
#(3)
#    model="DSLangmuir", param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#(3)                 model="BET", param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
    #, param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#    fill_value=df_gas1['p(mbar)'].max())

param_names = [param for param in gas1_isotherm.params.keys()]
print ('param_names:',param_names)
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
pyiast.plot_isotherm(gas2_isotherm)

gas2_isotherm.print_params()

#sys.exit()
gas2_interpolated = pyiast.InterpolatorIsotherm(
    df_gas2,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
    fill_value=df_gas2["loading(mmol/g)"].max())
# ## Perform IAST at same mixture conditions as binary GCMC simulations

# In[9]:
#iast_component_loadings = np.zeros((2, num))  # store component loadings here

nsteps=24

logpress=np.logspace(-3, 0, nsteps)
print (logpress)

for istep in range(nsteps):
        #total_pressure=0.1+istep*100
        total_pressure=10000*logpress[istep]
        print ("total pressure =", total_pressure)
	#mbar=1000*total_pressure
        iast_component_loadings = np.zeros((2, num))  # store component loadings here
        
        for i in range(num):
              y_gas2 = mol_frac[i]
              partial_pressures = total_pressure*np.array([y_gas2, 1.0 - y_gas2])
              iast_component_loadings[:, i] = pyiast.iast(
              partial_pressures, [gas2_isotherm, gas1_isotherm], verboseflag=False)
        
        fig = plt.figure()
        plt.plot(mol_frac, iast_component_loadings[0, :], color=color_key[gas2],
                       label='%s, IAST' %gas2) 
        plt.plot(mol_frac, iast_component_loadings[1, :], color=color_key[gas1],
                       label='%s, IAST' %gas1)

        plt.xlabel('Mole fraction %s in gas phase, $y_{%s}$' %(gas2,gas2))
        plt.ylabel('Gas uptake (loading(mmol/g))')

        plt.ylim(ymin=0)
        plt.xlim([0, 1.])

        plt.legend(loc='center left', prop={'size': 12})

        plt.tight_layout()
#plt.savefig("IAST_validation.pdf", format='pdf')
        plt.savefig("graphs/Graph_mixture_%s_%s_adsorption_%s_%d_mbar_at_%s.png" 
                               %(gas1, gas2, structure, total_pressure, temperature), 
                                                              format='png', dpi=250)
        plt.show()


        a['P=%s' %total_pressure]=iast_component_loadings[0,:]
        b['P=%s' %total_pressure]=iast_component_loadings[1,:]

#              bulk =(1-a['%s_mole_fract' %gas2])/a['%s_mole_fract' %gas2]
        selectivity['P=%s'%total_pressure]=b['P=%s' %total_pressure]/a['P=%s' %total_pressure]/bulk
#################################################################################################
file_name1="data/data_mixture_%s_adsorption_%s_at_%s.dat" %(gas2, structure, temperature)
file_name2="data/data_mixture_%s_adsorption_%s_at_%s.dat" %(gas1, structure, temperature)
file_name3="data/data_IAST_%s-%s_selectivity_%s_at_%s.dat" %(gas1,gas2,structure, temperature)
#a=pd.DataFrame(Loading1, columns=['#', 'ZIF69'])
a.to_csv(file_name1, sep='\t')
b.to_csv(file_name2, sep='\t')
selectivity.to_csv(file_name3, sep='\t')

# ## Another sanity check

# We use IAST for a three-component mixture of 5 bar methane (all the same!). This should yield the loading at 15 bar.

# In[11]:

iast_component_loadings = pyiast.iast(
    [5000., 5000., 5000.], [gas1_isotherm, gas1_isotherm, gas1_isotherm])
print("Sum of loadings: ", np.sum(iast_component_loadings))
print("Loading at 15 bar: ", gas1_isotherm.loading(15000.))
np.testing.assert_almost_equal(
    np.sum(iast_component_loadings), gas1_isotherm.loading(15000.))


# # Compare fits for ethane

# In[15]:

file_iso1="isotherms_out/out_%s_%s_isotherm_%s.dat" %(structure, gas1, temperature)
isotherm=pd.DataFrame(columns=['p(mbar)', 'loading(mmol/g)'])

gas1_interpolated = pyiast.InterpolatorIsotherm(
    df_gas1,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
    fill_value=df_gas1["loading(mmol/g)"].max())

p_plot = np.logspace(
    #np.log(df_gas1['p(mbar)'].min()), np.log(25000), num=200)
    np.log(0.1),np.log(150), num=200)

isotherm['p(mbar)']=p_plot
isotherm['loading(mmol/g)']=gas1_isotherm.loading(p_plot)
isotherm.to_csv(file_iso1, sep='\t', index=False)

fig = plt.figure()

plt.xlabel("Pressure (mbar)")
plt.ylabel("%s uptake (loading(mmol/g))" %gas1)

plt.scatter(
    df_gas1['p(mbar)'],
    df_gas1['loading(mmol/g)'],
    marker='o',
    color=color_key[gas1],
    s=40,
    label='Data',
    clip_on=False)

plt.plot(
    p_plot,
    gas1_interpolated.loading(p_plot),
    color='k',
    linewidth=1,
    label='Linear interpolation',
    linestyle='--')

plt.plot(
    p_plot,
    gas1_isotherm.loading(p_plot),
    color='k',
    linewidth=1,
    label='Langmuir fit')
plt.xscale("log")
plt.xlim(10, 2.5*10**4)
plt.ylim(0, df_gas1["loading(mmol/g)"].max()*1.5)

plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig("graphs/%s_isotherm_of_%s_at_%s.png" %(gas1, structure, temperature), format='png', dpi=250)

plt.show()


###################### second gas #############################
file_iso2="isotherms_out/out_%s_%s_isotherm_%s.dat" %(structure, gas2, temperature)
#isotherm=pd.DataFrame(columns=['p(mbar)' %gas2, 'loading(loading(mmol/g))' %gas2])

gas2_interpolated = pyiast.InterpolatorIsotherm(
    df_gas2,
    loading_key="loading(mmol/g)",
    pressure_key="p(mbar)",
    fill_value=df_gas2["loading(mmol/g)"].max())
    
#p_plot = np.logspace(
    #np.log(df_gas1['p(mbar)'].min()), np.log(25000), num=200)
#    np.log(7), np.log(15000), num=200)

isotherm['p(mbar)']=p_plot
isotherm['loading(mmol/g)']=gas2_isotherm.loading(p_plot)
isotherm.to_csv(file_iso2, sep='\t', index=False)

fig = plt.figure()

plt.xlabel("Pressure (mbar)")
plt.ylabel("%s uptake (loading(mmol/g))" %gas2)

plt.scatter(
    df_gas2['p(mbar)'],
    df_gas2['loading(mmol/g)'],
    marker='o',
    label='Data',
    color=color_key[gas2],
    s=40,
    clip_on=False)

plt.plot(
    p_plot,
    gas2_interpolated.loading(p_plot),
    color='k',
    linewidth=1,
    label='Linear interpolation',
    linestyle='--')

plt.plot(
    p_plot,
    gas2_isotherm.loading(p_plot),
    color='k',
    linewidth=1,
    label='Langmuir fit')
plt.xscale("log")
plt.xlim(10, 2.5*10**4)
plt.ylim(0, df_gas2["loading(mmol/g)"].max()*1.5)

plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig("graphs/%s_isotherm_of_%s_at_%s.png" %(gas2, structure, temperature), format='png', dpi=250)
plt.show()





# In[16]:

pyiast._MODELS
