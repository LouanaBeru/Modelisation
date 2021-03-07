##MODULES

import math as m
import numpy as np
import matplotlib.pyplot as plt

##DEFINITION VARIABLES
###DOMAINE
Lx = 1 #On se place sur [0,1]
Ly = 1 #On se place sur [0,1]

###MAILLAGE
Nx = 150 #en x
Ny = 150 #en y
x = np.linspace(0, Lx, Nx) #Nx points entre 0 et 1
y = np.linspace(0, Ly, Ny) #Ny points entre 0 et 1

###AUTRES
modOpti = 24 #mode optimal, trouve avec le projet 4 de modelisation
temps = 0.001 #Comme si on etait en t = 0

##CONSTRUCTION DE L INTEGRALE
###f_temps(x,y)    
def f(t, x, y):
    ft = m.exp( -(x-y)**2 / (4 * t) )
    return ft

###u(x,t)
cste = 1 / (2 * m.sqrt( m.pi * temps ) )

def u0(x):
    eval = m.exp( - 100 * (x - Lx / 2)**2 ) 
    return eval

###PRODUIT DE CONVOLUTION
u1 = []

for i in x:
    uy = 0
    for j in y:
        uy += f( temps, i, j ) * u0(j)
    u1.append( cste * uy )

plt.figure()
plt.plot(x, np.log(u1), color = 'red')

##AVEC np.convolve
###Vecteur f : évaluation de cste * f sur tous les points du maillage

fEval=[]

for i in x:
    fEval.append( cste * f(temps, i, 0) )

###Vecteur u0 : évaluation de u0 sur tous les points du maillage

u0Eval = []

for i in x:
    u0Eval.append( u0(i) )

###Produit de convolution

u2 = np.convolve( fEval, u0Eval, mode = 'full' )

###Graphe 

Nx2 = len(u2)
x2 = np.linspace(0, Lx, Nx2 )

plt.plot( x2, u2, color = 'blue')

##chaleur1dspec.py

#Initialisation 
k = 1
s= 100
Lx = 1
Nx = 150 #le maillage spatial sert a la representation de la solution et aux calculs des ps par integration numerique
hx = Lx/(Nx-1)
x = np.linspace(0,Lx,Nx)
modmax = modOpti
#introduire une boucle sur modmax et calculer l'erreur ||u(modmax)||
f=np.zeros(len(x))
u0=np.zeros(len(x))
u=np.zeros(len(x))
testx=np.zeros(modmax)
fbx=np.zeros((modmax,Nx))      #phi_n
fmp=np.zeros(modmax)           #projection sur phi_n de f
imp=np.zeros(modmax)           #projection sur phi_n de u0
#rhs
#Coeur du programme

for i in range(1,Nx):
    f[i]=0 #30*exp(-s*((x[i]-Lx/4)**2))  #(100*(x[i]**2)*((Lx-x[i])**2))/(Lx**4)
    u0[i]=m.exp(-s*((x[i]-Lx/2)**2)) #+exp(-2*s*((x[i]-Lx/3)**2))+exp(-3*s*((x[i]-2*Lx/3)**2))
    u[i]=u0[i]

#fb normalise
# 
for h in range(1,modmax):
    for i in range(1,Nx):
        fbx[h][i]=m.sin(m.pi*x[i]/Lx*h)               #phi_n
        testx[h]=testx[h]+fbx[h][i]*fbx[h][i]*hx  
    testx[h]=m.sqrt(testx[h])                       #norme L2 de phi_n
    for i in range(1,Nx):                         #normalisation
        fbx[h][i]=fbx[h][i]/testx[h]

        
#verifier l'orthonormalite des fbx ?
#<fbx[m],fbx[n]>=delta_mn


#projection f second membre et u0 condi init sur fbx

for h in range(1,modmax):
    for i in range(1,Nx):
        fmp[h]+=f[i]*fbx[h][i]*hx    # <f,\phi_n> = f_n 
        imp[h]+=u0[i]*fbx[h][i]*hx   # <u0,phi_n> = c_n

#somme serie

temps=0.001    #on doit retrouver la condition initiale
for i in range(1,Nx):
    u[i]=0
for h in range(1,modmax):
    al=(h**2)*(m.pi**2)/(Lx**2)*k
    coef=imp[h]*m.exp(-al*temps)
    for i in range(0,Nx-1):
         u[i]+=fbx[h][i]*coef
 
temps=0.02  #la solution a n'importe quel temps sans avoir a calculer les iter intermediaires
for i in range(1,Nx):
    u[i]=0
for h in range(1,modmax):
    al=(h**2)*(m.pi**2)/(Lx**2)*k
    coef=imp[h]*m.exp(-al*temps)
    coeff=fmp[h]*(1-m.exp(-al*temps))/al
    for i in range(0,Nx):
        u[i]+=fbx[h][i]*(coeff+coef)

plt.plot(x,u,'purple')

##GRAPHE FINAL 

plt.legend(('sans np.convolve', 'avec np.convolve', 'chaleur1dspec'))
plt.title('Convolution')
plt.show()


