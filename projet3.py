##README
#This code uses Gradient_Minimum, please make sure to save those two codes in the same directory before running projet3

###################################################################################################################################
##MODULES
import numpy as np
import matplotlib.pyplot as plt
from Gradient_Minimum import * #importe toutes les fonctions de ce programme

#u,t = k u,xx 
#u,t - k u,xx = f  equation de chaleur, mais ici, f =0

##CONDITIONS INITIALES
###PHYSICAL PARAMETERS
Rapport_temp = 4 #ratio temp entre x=0 et L, les quantites sont sans dimension (voir chap. 3)
K = 0.5    #Reference Diffusion coefficient
L = 1.0     #Domain size
Time = 1  #Integration timeRapport_temp=4 #ratio temp entre x=0 et L, les quantites sont sans dimension (voir chap. 3)

###NUMERICAL PARAMETERS
NX = 100  #Number of grid points
ifre=10
eps=0.01
dx = L/(NX-1)  #Grid step (space)
dt=dx**2/(2*K) #CFL stability condition
NT = int(Time/dt)   #Number of time steps necessary to reach Time under stability condition
CNTRL = 10

##MAIN PROGRAM
###Initialisation
x = np.linspace(0.0,1.0,NX) #La porte fait 10 cm d'épaisseur, en 0, on est à l'intérieur, en 1, on est à l'extérieur ( du côté du feu)
T = np.zeros((NX)) #np.sin(2*np.pi*x) #Température en fonction de là où on se place dans l'espace
T0 = 1  
T[0] = T0 #Température à l'intérieur
T[NX-1] = T0*Rapport_temp #Température à l'extérieur = là où il y a le feu, il fait Rapport_temp fois plus chaud

rest = []
RHS = np.zeros((NX))

###Main loop en temps
n=0
res=1
res0=1
for n in range(0,NT):#On doit avoir fini avant NT

    while(n<NT and res/res0>eps):
        n+=1
###discretization of the second order derivative (Laplacian)
        res=0
        for j in range (1, NX-1):
            RHS[j] = dt*(K*(T[j-1]-2*T[j]+T[j+1])/(dx**2)+1)
            res+=abs(RHS[j])

        for j in range (1, NX-1):
            T[j] += RHS[j]
            RHS[j]=0


        if (n == 1 ):
            res0=res

        rest.append(res)
    
###Ajout de l'isolant
#Fonction : KX
#Parametres : x : list / maillage, CNTRL : float / plus il est grand, plus c'est isole
#Retourne : float / la temperature avec l'isolant
def KX(x, CNTRL):
    temp = K /(1+CNTRL*np.exp(-1.e6*(x-0.5)**8)) #température avec l'isolant
    return temp

#Fonction : KXFixCNTRL
#Parametres : CNTRL : float / plus il est grand, plus c'est isole
#Retourne : la fonction KX ou CNTRL est fixe ( il vaut la valeur qu'on donne a la fonction)
def KXFixCNTRL(CNTRL):
    return lambda x: KX(x, CNTRL)

#Fonction : KXFixX
#Parametres : x : list / maillage
#Retourne :  la fonction KX ou x est fixe ( il vaut la valeur qu'on donne a la fonction)
def KXFixX(x):
    return lambda CNTRL: KX(x, CNTRL)

#Fonction : J
#Parametres : CNTRL : float / plus il est grand, plus c'est isole
#Retourne : cost : float / le cout
def J(CNTRL):
    gradKX = Gradient(KXFixX(0.5), 0, L, NX, eps)
    cost = abs(gradKX[CNTRL] - 0.5*gradKX[0])
    return cost

#Fonction : cost
#Parametres : /
#Retourne : une liste du cout pour differentes valeurs de CNTRL
def cost():
    cost = []
    for i in range(NX):
        CNTRL = i
        cost += [J(i)]
    return cost

##AFFICHAGE
plt.figure()      
plt.plot(x, cost(), color = 'orange')

kx = np.zeros(len(x))
for i in range(len(x)):
    kx[i] = K /(1+4.141*np.exp(-1.e6*(x[i]-0.5)**8))
plt.figure()
plt.plot(x,np.log10(kx), color = 'blue' )
plt.plot(x, T , color = 'yellow')
plt.show()

#CNTRL divise par deux la propagation de la chaleur -> comment le représenter ?
