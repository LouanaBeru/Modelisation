import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import cm
from array import array


def tri(list):
    list.sort()
    return list

def norm2(list):
    norm2 = 0
    for i in range(len(list)):
        norm2 = norm2 + list[i]**2
    return np.sqrt(norm2)

##VARIABLES
L1=[]
L2=[]
N=[]
NL1=[]
NL2=[]

# génère un entier au hasard entre 3 et 200 et on met les différentes valeurs prises dans une liste
N=random.choices(range(3, 201), k=20) 
N=tri(N) # on tri par ordre croissant ( pour l'affichage)


for i in range(20) : #on le fait 20 fois
    # NUMERICAL PARAMETERS
    NX = N[i] #ieme nombre de points du maillage ( que je noterai NXi dans les légendes)
    L=2*3.141592   #2*math.pi
    dx = L/(NX-1)  #Grid step (space)

    # Initialisation
    x = np.linspace(0.0,L,NX) #x = [x_0, ..., x_NXi] -> la liste des points du maillage

    T = np.sin(x)
    Txe = np.cos(x)
    Txxe = -np.sin(x)

    Tx = np.zeros((NX))
    Txx = np.zeros((NX))

    #discretization of the second order derivative (Laplacian)
    for j in range (1, NX-1):
        Tx[j] = (T[j+1]-T[j-1])/(dx*2)
        Txx[j] = (T[j-1]-2*T[j]+T[j+1])/(dx**2)
        
    #Tx and Txx on boundaries
    # use extrapolation in order to have (T,x),xx=0
    #(T,x),xx= (Tx0 -2*Tx1+Tx2) =0
    Tx[0] = 2*Tx[1]-Tx[2]
    Tx[NX-1] = 2*Tx[NX-2]-Tx[NX-3]
    Txx[0] = 2*Txx[1]-Txx[2]
    Txx[NX-1] = 2*Txx[NX-2]-Txx[NX-3]

    L1= abs(Tx-Txe) #liste de liste des erreurs en chaque point du maillage pour la première dérivée
    L2= abs(Txx-Txxe) #liste de liste des erreurs pour la deuxième dérivée
    NL1 = NL1 + [norm2(L1)]
    NL2 = NL2 + [norm2(L2)]

plt.figure(1) #pour afficher
plt.xlabel(u'$NX$', fontsize=26)
plt.ylabel(u'$Erreur$', fontsize=26, rotation=90)
plt.plot(N,np.log10(NL1),label='L1',color=(1,0,0,1)) #On trace l'erreur (abs(Tx-Txe) pour tous les points xi pour la couleur, il fallait que la couleur soit moins intense quand il y a moins de points ( la couleur s'exprime en RGBA, A étant la transparence) et c'est un pourcentage ( donc un nbr entre 0 et 1)
plt.title(u'Evolution de l\' erreur pour la première dérivée')   
    
plt.figure(2) #pour afficher
plt.xlabel(u'$NX$', fontsize=26)
plt.ylabel(u'$Erreur$', fontsize=26, rotation=90)
plt.plot(N,np.log10(NL2),label='L2',color=(0,0,1,1)) 
plt.title(u'Evolution de l\' erreur pour la deuxième dérivée')
plt.show()


#ameliorer ce code en ajoutant les fonctionalites suivantes:
# parameterer NX et etudier l'evolution de l'erreur L1 pour 20 valeurs de NX de 3 a 200
# tracer les courbes de l'erreur en fonction de NX en echelle log10


