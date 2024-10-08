# coding: utf-8

# # Test

# We test the IAST code with a binary mixture of CO2 and N2 (or O2 and N2) adsorbed in four different samples(struct-list).
#
# The isotherms to be read ar the 16 interpolated isotherms at temperatures from 278 to 313K.
# Interpoltated pure-component adsorption isotherms from 278 to 313K K are present in the files:
# * isosteric/interpolated_%s_%s_isotherms_on_%s.dat
# the selectivities (at variable temperature and pressure are stored in data/data_IAST_%s-%s_selectivity_%s_at_%s.dat
# 
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
gas1=gas_list[2]
gas2=gas_list[1]

structure=struct_list[3]
temperature=temp_list[0]

press_list=[100,1000,10000]

total_pressure=press_list[2]


print(gas1,gas2)
# ## Load pure-component isotherm data as Pandas dataframes

# In[2]:
x0=283
xn=313
nT=16
Textra=np.linspace(x0,xn,nT)
num_plots=nT



#
file1="isosteric/interpolated_%s_%s_isotherms_on_%s.dat" %(num_plots, gas1, structure) 

df_gas1 = pd.read_csv(file1, sep='\s+')
df_gas1.head()
icol=10  #a random temperature step
temperature=Textra[icol]
print (df_gas1)
print (df_gas1[str(Textra[icol])])

# ### Ethane

# In[3]:

file2="isosteric/interpolated_%s_%s_isotherms_on_%s.dat" %(num_plots, gas2, structure)
df_gas2 = pd.read_csv(file2, sep='\s+')
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

#fit model parameters record for isotherms 1 and 2 at variable temp
column_names = [Textra[icol] for icol in range(len(Textra))]
parameters1=pd.DataFrame(index=['M','Ka','Kb'], columns=column_names)
parameters2=pd.DataFrame(index=['M','K'], columns=column_names)


#partial_pressures=mol_fract*total_pressure

print(total_pressure)

for icol in range(len(Textra)):
    fig, ax = plt.subplots()

    plt.scatter(
#    df_gas2['p(mbar)'],
#    df_gas2['loading(mmol/g)'],
          df_gas2[str(Textra[icol])],
          df_gas2['mmol/g'],
          label=gas2,
          color=color_key[gas2],
          s=20,
          marker='s')
    
    plt.scatter(
          df_gas1[str(Textra[icol])],
          df_gas1['mmol/g'],
          label=gas1,
          color=color_key[gas1],
          s=20,
          marker='o')

    plt.xlabel('Pressure (mbar)')
    plt.ylabel('Gas uptake (loading(mmol/g))')

    plt.xlim([0, 20000])
    plt.ylim(ymin=0)

    plt.legend(loc='lower right')

    plt.tight_layout()
   #plt.savefig("pure_component_isotherms.pdf", format='pdf')
#    plt.savefig("iso_plots/pure_component_%s_and_%s_isotherms_of_%s_at_%s.png"  
#           %(gas1, gas2, structure, str(Textra[icol])), format='png', dpi=250)
    plt.show()
###############
#sys.exit()
# ### Interpolator isotherm for Methane

# In[7]:
    model1="DSLangmuir"
#    model1="BET"
    model1="Langmuir"
    
    
    gas1_isotherm = pyiast.ModelIsotherm(
#gas1_isotherm = pyiast.InterpolatorIsotherm(
         df_gas1,
         loading_key="mmol/g",
         pressure_key=str(Textra[icol]),
#    fill_value=df_gas1["loading(mmol/g)"].max())
#    model="Langmuir")
#    model="TemkinApprox", param_guess= {'M':10, 'K':0.03,'theta':0.0001})
#(0) 
#         model=model1, param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
#(1) (default)      
    model=model1)
    #, param_guess= {'M':2.1, 'K':0.001180,'Kb':0.000024})
#    model="BET")
#(2)             
#    model=model1, param_guess= {'M1':3.36, 'K1':0.002,'M2':3.2,'K2':1.78})
#(3)
#    model="DSLangmuir", param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#(3)                 model="BET", param_guess= {'M':2.1, 'Ka':0.001180,'Kb':0.000024})
    #, param_guess= {'M1':100, 'K1':0.1,'M2':100,'K2':0.1})
#    fill_value=df_gas1['p(mbar)'].max())
    param_names = [param for param in gas1_isotherm.params.keys()]
    print ('param_names for gas1:',param_names)
# ### Interpolator isotherm for ethane
    gas1_isotherm.print_params()
    
    if model1 == "Langmuir":
       parameters1.at['M',Textra[icol]]=gas1_isotherm.params['M']
       parameters1.at['K',Textra[icol]]=gas1_isotherm.params['K']
    if model1 == "BET":
       parameters1.at['M',Textra[icol]]=gas1_isotherm.params['M']
       parameters1.at['Ka',Textra[icol]]=gas1_isotherm.params['Ka']
       parameters1.at['Kb',Textra[icol]]=gas1_isotherm.params['Kb']
    if model1 == "DSLangmuir":
       parameters1.at['M1',Textra[icol]]=gas1_isotherm.params['M1']
       parameters1.at['K1',Textra[icol]]=gas1_isotherm.params['K1']
       parameters1.at['M2',Textra[icol]]=gas1_isotherm.params['M2']
       parameters1.at['K2',Textra[icol]]=gas1_isotherm.params['K2']
#    pyiast.plot_isotherm(gas1_isotherm)
# In[8]:
####******************************************************************************************************66
#sys.exit()
####******************************************************************************************************77
    model2="Langmuir"
    gas2_isotherm = pyiast.ModelIsotherm(
#gas2_isotherm = pyiast.InterpolatorIsotherm(
           df_gas2,
           loading_key="mmol/g",
           pressure_key=str(Textra[icol]),
#    fill_value=df_gas2["loading(mmol/g)"].max())
#    model="Langmuir")
           model=model2)
    #model="DSLangmuir", param_guess= {'M1':1300, 'K1':10,'M2':1000,'K2':10})
#    fill_value=df_gas2["p(mbar)"].max())
    param_names = [param for param in gas2_isotherm.params.keys()]
#    pyiast.plot_isotherm(gas2_isotherm)
    print ('param_names for gas2:',param_names)
    gas2_isotherm.print_params()
    
#    parameters2=pd.DataFrame(index=['M','K'],columns=Textra.keys())
    parameters2.at['M',Textra[icol]]=gas2_isotherm.params['M']
    parameters2.at['K',Textra[icol]]=gas2_isotherm.params['K']

# In[9]:
#iast_component_loadings = np.zeros((2, num))  # store component loadings here
    isotherm1=pd.DataFrame(columns=['p(mbar)', 'loading(mmol/g)'])
    isotherm2=pd.DataFrame(columns=['p(mbar)', 'loading(mmol/g)'])
    p_plot = np.logspace(np.log(0.1),np.log(150), num=200)
    
    isotherm1['p(mbar)']=p_plot
    isotherm1['loading(mmol/g)']=gas1_isotherm.loading(p_plot)
    isotherm2['p(mbar)']=p_plot
    isotherm2['loading(mmol/g)']=gas2_isotherm.loading(p_plot)
#isotherm.to_csv(file_iso1, sep='\t', index=False)

    fig = plt.figure()

    plt.xlabel("Pressure (mbar)")
    plt.ylabel("%s uptake (loading(mmol/g) at T = %s K)" %(gas1, str(Textra[icol])))

    plt.scatter(
        df_gas1[str(Textra[icol])],
        df_gas1['mmol/g'],
        marker='o',
        color=color_key[gas1],
    #s=40,
        label=gas1,
        clip_on=False)
    
    plt.scatter(
        df_gas2[str(Textra[icol])],
        df_gas2['mmol/g'],
        marker='o',
        color=color_key[gas2],
    #s=40,
        label=gas2,
        clip_on=False)


    plt.plot(
        p_plot,
        gas1_isotherm.loading(p_plot),
        color='k',
        linewidth=1,
        label='models fit' )
    plt.plot(
        p_plot,
        gas2_isotherm.loading(p_plot),
        color='k',
        linewidth=1)
    
    plt.xscale("log")
    plt.xlim(10, 2.5*10**4)
    plt.ylim(0, df_gas1["mmol/g"].max()*1.5)

    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig("iso_plots/graph_%s_and_%s_isotherms_of_%s_at_%s.png" %(gas1, gas2, structure, str(Textra[icol])), format='png', dpi=250)

    plt.show()
####################################################################################################################################
###############              The single component isotherms at the specified temperature are set          ##########################
###############               now we start our IAST analysis at fixed pressure and varying T              ##########################
####################################################################################################################################    
#    total_pressure=10000*logpress[istep]
    print ("total pressure =", total_pressure)
############                   Initialize iast_component_loadings                      #############################################    
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
    plt.ylabel("gas uptake (mmol/g) at T = %s K)" %str(Textra[icol]))

    plt.ylim(ymin=0)
    plt.xlim([0, 1.])

    plt.legend(loc='center left', prop={'size': 12})

    plt.tight_layout()
#plt.savefig("IAST_validation.pdf", format='pdf')
    plt.savefig("iso_plots/Graph_mixture_%s_%s_adsorption_%s_%d_mbar_at_%s.png" 
                               %(gas1, gas2, structure, total_pressure, str(Textra[icol])), 
                                                              format='png', dpi=250)
    plt.show()          
    a['T=%s' %str(Textra[icol])]=iast_component_loadings[0,:]
    b['T=%s' %str(Textra[icol])]=iast_component_loadings[1,:]              
    selectivity['T=%s' %str(Textra[icol])]=(iast_component_loadings[1,:]/iast_component_loadings[0,:])/bulk

###########                     pressure is in mbar #########################################

paramFile1="iso_params/params_%s_%s_%s.dat" %(model1, gas1, structure)
parameters1.to_csv(paramFile1,sep='\t')
paramFile2="iso_params/params_%s_%s_%s.dat" %(model2, gas2, structure)
parameters2.to_csv(paramFile2,sep='\t')

file_name1="variableT_out/data_mixture_%s_adsorption_%s_at_%s_mbar.dat" %(gas2, structure, total_pressure)
file_name2="variableT_out/data_mixture_%s_adsorption_%s_at_%s_mbar.dat" %(gas1, structure, total_pressure)
file_name3="variableT_out/data_IAST_%s-%s_selectivity_%s_at_%s_mbar.dat" %(gas1,gas2,structure, total_pressure)
#a=pd.DataFrame(Loading1, columns=['#', 'ZIF69'])
a.to_csv(file_name1, sep='\t')
b.to_csv(file_name2, sep='\t')
selectivity.to_csv(file_name3, sep='\t')
sys.exit()


# In[16]:

pyiast._MODELS
