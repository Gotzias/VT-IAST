import pyiast
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import exp
import sys
from scipy import stats
from scipy.interpolate import interp1d

# Matplotlib settings
plt.style.use('bmh')


###############################################################################################################################
###############################################################################################################################
#######################        This code reads the single component adsorption isotherm      ##################################
######################    of either N2, CO2 or O2  at thee temperatures such as 283, 298 and 313   ############################
#######################  It computes the isosteric heat (DH) of adsorption of the adsorbate gas    ############################
########################  Then it creates a series of 16 interpolated isotherms             ###################################
########################                   ranging from 283K to 313 K                        ##################################
########################     The isotherms to be read are stored in folder heat/             ##################################
###############################################################################################################################
###############################################################################################################################
color_key = {'T1':'b', 'T2':'r', 'T3':'g', 'T4':'orange', 'T5':'maroon', 'T6':'grey'}
color_key = {'all':'b', 'kno':'r', 'gme':'g'}
ArrayofGases=['CH4', 'CO2', 'H2']
ArrayofStructs=['MgBH', 'ZIF72', 'mgf', 'cli']
Arrayofranges=['section_1', 'section_2', 'section_3', 'section_4', 'section_5']
#ArrayofTemp=['078', '083', '088', '093', '098', '100', '125', '150', '175', '190', '250', '290']
#ArrayofTemp=['078', '083', '088', '093', '098', '100']
#ArrayofTemp=['100', '125', '150', '175', '190']
#ArrayofTemp=['190','250']
#ArrayofTemp=['250','290']

color_key = {'CO2': 'g', 'N2': 'b', 'O2': 'r'}
struct_list = ['Alumina', 'WS2050', 'XLBFK-12', 'XLBFK-16']
temp_list=['283K', '298K', '313K']
temp2_num=['283', '298', '313']
gas_list=['CO2', 'N2', 'O2']
gas1=gas_list[1]


structure=struct_list[0]
temperature=temp_list[0]

Gas=ArrayofGases[2]
Str=ArrayofStructs[0]
rangeis=Arrayofranges[0]

modelname='Langmuir'
#modelname='BET'
#modelname='DSLangmuir'
print("reading the %s isotherms on %s with %s-fit" %(gas1, structure, modelname))
df_data=pd.read_csv("heat/isos_%s_%s_%s.dat" %(modelname, gas1, structure), delimiter='\t')
df_data.head()

ArrayofTemps=list(df_data)
Tinv=[float(i) for i in ArrayofTemps[1:]]
TK=pd.Series(Tinv)
Tinv= np.reciprocal(Tinv)
print('****')
print(Tinv)
#sys.exit()
P_500=df_data['mbar']
data_give=pd.DataFrame(P_500, columns=['mbar'])
print(data_give['mbar'])
#sys.exit()
#P_2000=df_data2['kPa']
#data_give2=pd.DataFrame(P_2000, columns=['kPa'])
#Find the range of uptakes for the qst plots (i.e q=[max(qo),min(qn)]) 

qn=10000
q0=-100
df_nonzero = df_data[df_data != 0.0]
for icol in range(1,len(df_data.columns)):
    q=df_data[ArrayofTemps[icol]].max()  
    qq=df_nonzero[ArrayofTemps[icol]].min()
    print(qq,q)
    if q<qn : qn=q
    if q0<qq:q0=qq

#see q=[max(qo),min(qn)])
print('q=[%s,%s]' %(q0,qn))

y=df_data['mbar']
q_heat=np.linspace(q0,qn,1000)
Matrix=pd.DataFrame(q_heat, columns=['mmol/g'])
Coeff=pd.DataFrame(q_heat, columns=['mmol/g'])
Interp_isos=pd.DataFrame(q_heat, columns=['mmol/g'])

for icol in range(1,len(df_data.columns)):
    x=df_data[ArrayofTemps[icol]]
    f=interp1d(x,y)
    Matrix[ArrayofTemps[icol]]=np.log(f(q_heat))

temp_frame=[]
temp_frame2=[]
#x=[Tinv[0],Tinv[1],Tinv[2], Tinv[3], Tinv[4],Tinv[5]]
x=[Tinv[0],Tinv[1],Tinv[2]]
#sys.exit()
print("Tinv:",Tinv)
print("x",x)

#x=[Tinv[0],Tinv[1]]
for irow, row in Matrix.iterrows():
#    y=[row[ArrayofTemps[1]],row[ArrayofTemps[2]],row[ArrayofTemps[3]],row[ArrayofTemps[4]],row[ArrayofTemps[5]],row[ArrayofTemps[6]]]
    y=[row[ArrayofTemps[1]],row[ArrayofTemps[2]],row[ArrayofTemps[3]]]
    slope,intercept, r_value, p_value, std_err=stats.linregress(x,y)
    temp_frame.append(slope)
    temp_frame2.append(intercept)
    print(irow,":",y,"slope",slope,"intercept",intercept)

#extrapolate isothemrs for many Temperatures (Textra) in the form q,p 
#for q[q0,qmax] and lnp=dh/(RT)+c
#sys.exit()
x0=283
xn=313
nT=16
Textra=np.linspace(x0,xn,nT)
Textra_inv=1./Textra
Coeff['slope']=temp_frame
Coeff['Intrcpt']=temp_frame2

Previous_coef_file="heat/params_%s_%s_%s.dat" %(modelname, gas1, structure)
quad_coef=pd.read_csv(Previous_coef_file, delimiter='\t')
quad_coef.head()
#sys.exit()
if modelname == "Langmuir":
   mu1=[quad_coef[i][0] for i in temp_list]
   kappa_a=[quad_coef[i][1] for i in temp_list]
if modelname == "BET":
   mu1=[quad_coef[i][0] for i in temp_list]
   kappa_a=[quad_coef[i][1] for i in temp_list]
   kappa_b=[quad_coef[i][2] for i in temp_list]
if modelname == "DSLangmuir":
   kappa_a=[quad_coef[i][0] for i in temp_list]
   kappa_b=[quad_coef[i][1] for i in temp_list]
   mu1=[quad_coef[i][2] for i in temp_list]
   mu2=[quad_coef[i][3] for i in temp_list]
   
#kappas=[quad_coef[i][0] for i in temp_list]
#saturation_load=[quad_coef[i][1] for i in temp_list]
#alphas=[quad_coef[i][2] for i in temp_list]
#betas=[quad_coef[i][3] for i in temp_list]
#ciis=[quad_coef[i][4] for i in temp_list]
x_param=[float(i) for i in ArrayofTemps[1:]]
print ("xparam=", x_param)
print ("it's %s", modelname)  
if modelname == "Langmuir": print ("saturation_loading=", mu1, "K = ", kappa_a)
if modelname == "BET": print ("saturation_loading=", mu1, "K1=", kappa_a,"K2=", kappa_b)
if modelname == "DSLangmuir": print ("saturation_loading-a = ", mu1, "K1 = ", kappa_a,"saturation_loading-b = ", mu2, "K2 = ", kappa_b)

x_new=[float(i) for i in Textra]
#print("kappas= ", kappas)
#print("x_param= ", x_param)
if modelname == "Langmuir":
   f_kappa=interp1d(x_param,kappa_a,kind='quadratic')
   f_sat_load=interp1d(x_param,mu1,kind='quadratic')
if modelname == "BET":
   f_kappa=interp1d(x_param,kappa_a,kind='quadratic')
   f_sat_load=interp1d(x_param,mu1,kind='quadratic')
   f_kappa2=interp1d(x_param,kappa_b,kind='quadratic' )
if modelname == "DSLangmuir":
   f_kappa=interp1d(x_param,kappa_a,kind='quadratic')
   f_sat_load=interp1d(x_param,mu1,kind='quadratic')
   f_kappa2=interp1d(x_param,kappa_b,kind='quadratic' )
   f_sat_load2=interp1d(x_param,mu2,kind='quadratic' )

if modelname == "Langmuir":
   plt.plot(x_param, kappa_a, 'o', x_new, f_kappa(x_new), '-')
   plt.legend(['data', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, mu1, 'o', x_new, f_sat_load(x_new), '-')
   plt.legend(['data', 'quadratic'], loc='best')
   plt.show()

if modelname == "BET":
   plt.plot(x_param, kappa_a, 'o', x_new, f_kappa(x_new), '-')
   plt.legend(['K_BET_a', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, mu1, 'o', x_new, f_sat_load(x_new), '-')
   plt.legend(['Saturation loading', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, kappa_b, 'o', x_new, f_kappa2(x_new), '-')
   plt.legend(['K_BET_b', 'quadratic'], loc='best')
   plt.show()

if modelname == "DSLangmuir":
   plt.plot(x_param, kappa_a, 'o', x_new, f_kappa(x_new), '-')
   plt.legend(['K_1', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, mu1, 'o', x_new, f_sat_load(x_new), '-')
   plt.legend(['M_1', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, kappa_b, 'o', x_new, f_kappa2(x_new), '-')
   plt.legend(['K_b', 'quadratic'], loc='best')
   plt.show()
   plt.plot(x_param, mu2, 'o', x_new, f_sat_load2(x_new), '-')
   plt.legend(['M_2', 'quadratic'], loc='best')
   plt.show()
#plt.plot(x_param, ciis, 'o', x_new, f_cii(x_new), '-')
#plt.legend(['data', 'cubic'], loc='best')
#plt.show()


num_plots=nT
#colormap=plt.cm.gist_ncar
#plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
labels=[]
for icol in range(len(Textra)):
    Interp_isos[Textra[icol]]=Coeff['slope']*Textra_inv[icol]+Coeff['Intrcpt']
    Interp_isos[Textra[icol]]=np.exp(Interp_isos[Textra[icol]])
    plt.plot(Interp_isos[Textra[icol]],Interp_isos['mmol/g'])
    labels.append('T=%dK' %Textra[icol])


plt.legend(labels, ncol=4, loc='lower center', 
           bbox_to_anchor=[0.5, 0.0], 
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)
plt.savefig("%s_%s_extrapolation.png" %(Gas, Str), format='png', dpi=250)
plt.show()
 
temp_new=[-8.314*x/1000. for x in temp_frame]
Matrix['slope']=temp_frame
Matrix['DH']=temp_new
fileall="isosteric/experimental_%s_isotherms_on_%s.dat" %(gas1, structure) 
Matrix.to_csv(fileall, index= False, sep='\t')
file2="isosteric/interpolated_%s_%s_isotherms_on_%s.dat" %(num_plots, gas1, structure)    
Interp_isos.to_csv(file2, index=False, sep='\t')   

