import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import pandas as pd
from io import StringIO
import random
import scipy
from scipy import stats
import statsmodels.api as sm


golden_mean = (math.sqrt(5)-1.0)/2.0       # Aesthetic ratio
fig_width = 7+3/8                          # width  in inches
fig_height = fig_width*golden_mean         # height in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
          'axes.titlesize': 18,
          'axes.labelsize': 19,
          'axes.linewidth': 0.7, 
          'axes.grid': False,
          'axes.labelweight': 'normal',  
          'font.family': 'serif',
          'font.size': 18.0,
          'font.weight': 'normal',
          'text.color': 'black',
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': True,
          'legend.fontsize': 16,
          'figure.dpi': 700,
          'figure.figsize': fig_size,
          'savefig.dpi': 700,
         }

pylab.rcParams.update(params)

#################### IMPORTACION DE DATA ####################

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G_Y2H = nx.read_edgelist('yeast_Y2H.txt')
G_LIT = nx.read_edgelist('yeast_LIT.txt')

G_LIT_Reguly = nx.Graph()
data_lit_reguly = np.genfromtxt("yeast_LIT_Reguly.txt",dtype='unicode',delimiter = '\t')
yeast_LIT_Reguly_raw = [data_lit_reguly[1:,0],data_lit_reguly[1:,1]]

i=1
while i<len(data_lit_reguly[1:,0]):
    G_LIT_Reguly.add_edge(data_lit_reguly[i,0],data_lit_reguly[i,1])
    i+=1


G_ESSENTIAL = nx.Graph()
data_G_ESSENTIAL = np.genfromtxt("Essential_ORFs_paperHe.txt",dtype='unicode',delimiter = '\t')
yeast_G_ESSENTIAL_raw = [data_G_ESSENTIAL[2:,1]]

i = 0
l = len(yeast_G_ESSENTIAL_raw[0][:])

while i < l:
    proteina = yeast_G_ESSENTIAL_raw[0][i]
    yeast_G_ESSENTIAL_raw[0][i] = proteina.rstrip(" ")
    i+=1


#################### Agrega esencialidad ####################
def agrega_esencialidad(red,lista_esencial):
    
    j = 0
    l = len(red)
    while j < l:
        nombre = list(red)[j]
        red.node[nombre]['esencialidad']=0
        j += 1
    
    for i in lista_esencial[0]:
        if (i in red):
            red.node[i]['esencialidad']=1
    return

agrega_esencialidad(G_APMS,yeast_G_ESSENTIAL_raw)
agrega_esencialidad(G_Y2H,yeast_G_ESSENTIAL_raw)
agrega_esencialidad(G_LIT,yeast_G_ESSENTIAL_raw)
agrega_esencialidad(G_LIT_Reguly,yeast_G_ESSENTIAL_raw)


def maximo_grado(red):
    cantidad_nodos = red.number_of_nodes()
    maximo = 0
    i=0
    while i<cantidad_nodos:
        if list(red.degree)[i][1]>maximo:
            maximo = list(red.degree)[i][1]
        i+=1
    return maximo

def he(red):

	k_max = maximo_grado(red)
	vector_esencialidad_nodo = np.zeros(k_max+1)
	histo_k = np.zeros(k_max+1)
	vector_k = np.linspace(0,k_max,k_max+1)
	#print(vector_k)
	pe = np.zeros(k_max+1)
	y_axis = np.zeros(k_max+1)

	n = len(red)
	i = 0
	while i < n:
		k = list(red.degree)[i][1]
		histo_k[k] += 1
        
		nombre = list(red)[i]
		esencialidad_nodo = red.node[nombre]['esencialidad']
        
		if (esencialidad_nodo==1):
			vector_esencialidad_nodo[k] += 1
		i+=1

	for i in range (0,len(vector_esencialidad_nodo)):
		if histo_k[i]>0:
			pe[i] = vector_esencialidad_nodo[i]/histo_k[i]
			if(i<11):
				y_axis[i] = math.log(1-pe[i])
	return (vector_k,pe,y_axis)

#print(he(G_LIT))

resultado_he_LIT = he(G_LIT)
resultado_he_LIT_Reguly = he(G_LIT_Reguly)
resultado_he_APMS = he(G_APMS)
resultado_he_Y2H = he(G_Y2H)

t = np.linspace(0,10,100)

'''
coef_fit_LIT = np.polyfit(resultado_he_LIT[0][1:10],resultado_he_LIT[2][1:10],1)
fit_LIT = coef_fit_LIT[1]+t*coef_fit_LIT[0]
coef_fit_LIT_Reguly = np.polyfit(resultado_he_LIT_Reguly[0][1:10],resultado_he_LIT_Reguly[2][1:10],1)
fit_LIT_Reguly = coef_fit_LIT_Reguly[1]+t*coef_fit_LIT_Reguly[0]


plt.plot(resultado_he_LIT[0][1:10],resultado_he_LIT[2][1:10],'^g',label='LIT')
plt.plot(t,fit_LIT,'-g')
plt.plot(resultado_he_LIT_Reguly[0][1:10],resultado_he_LIT_Reguly[2][1:10],'^y',label='LIT Reguly')
plt.plot(t,fit_LIT_Reguly,'-y')


coef_fit_APMS = np.polyfit(resultado_he_APMS[0][1:10],resultado_he_APMS[2][1:10],1)
fit_APMS = coef_fit_APMS[1]+t*coef_fit_APMS[0]
coef_fit_Y2H = np.polyfit(resultado_he_Y2H[0][1:10],resultado_he_Y2H[2][1:10],1)
fit_Y2H = coef_fit_Y2H[1]+t*coef_fit_Y2H[0]

plt.plot(resultado_he_APMS[0][1:10],resultado_he_APMS[2][1:10],'^r',label='AP/MS')
plt.plot(t,fit_APMS,'-r')
plt.plot(resultado_he_Y2H[0][1:10],resultado_he_Y2H[2][1:10],'^k',label='Y2H')
plt.plot(t,fit_Y2H,'-k')



pylab.xlabel('protein connectivity (k)')
pylab.ylabel('ln(1-P$_{E}$)')
pylab.ylim(-1.2, 0)
pylab.xlim(0, 11)
plt.legend()
pylab.savefig('fig2_he_b2.eps', format='eps', dpi=300, bbox_inches='tight')

'''
'''
plt.plot(resultado_he[0][1:10],resultado_he[1][1:10],'og',label='LIT Reguly')
pylab.xlabel('protein connectivity (k)')
pylab.ylabel('Pe')
#pylab.ylim(0, 1)
#pylab.xlim(0, 1)
plt.legend()
pylab.savefig('fig2_he_a.eps', format='eps', dpi=300, bbox_inches='tight')


##### Resultados Fit #########

slope, intercept, r_value, p_value, std_err = stats.linregress(resultado_he_APMS[0][1:10],resultado_he_APMS[2][1:10])
print "APMS",round(slope,3), round(intercept,3), round(r_value*r_value,3)

slope, intercept, r_value, p_value, std_err = stats.linregress(resultado_he_LIT[0][1:10],resultado_he_LIT[2][1:10])
print "LIT",round(slope,3), round(intercept,3), round(r_value*r_value,3)

slope, intercept, r_value, p_value, std_err = stats.linregress(resultado_he_LIT_Reguly[0][1:10],resultado_he_LIT_Reguly[2][1:10])
print "LIT Reguly",round(slope,3), round(intercept,3), round(r_value*r_value,3)

slope, intercept, r_value, p_value, std_err = stats.linregress(resultado_he_Y2H[0][1:10],resultado_he_Y2H[2][1:10])
print "Y2H",round(slope,3), round(intercept,3), round(r_value*r_value,3)

'''
##### Resultados segundo Fit #########

## Sacado de https://stackoverflow.com/questions/26050855/getting-uncertainty-values-in-linear-regression-with-python


x = resultado_he_Y2H[0][1:10]
y = resultado_he_Y2H[2][1:10]


x = sm.add_constant(x)

model = sm.OLS(y,x)
results = model.fit()


print(results.params)
print(results.bse )
