# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 08:57:18 2021

@author: Chapati
"""

#numpy.gradient(f,varargs, axis, edge_order) -> donne le gradient d'une fonction

#Se placer sur un intervalle et définir un nombre de points, créer le vecteur x, qui a pour coordonnées les différents points
# ensuite, demander à l'utilisateur de rentrer une fonction f en fonction de x
#créer un vecteur f qui a pour coordonnées f appliquée aux x[i]
# => on aura les valeurs de f sur tous nos points en discret


# calculer le gradient, div, laplacien par le meme principe

# affichage : relier les points
import numpy as np
import matplotlib.pyplot as plt
import random

L = 5 #longueur de l'axe X, on est sur ]0;5[
NX = 100 #Nombre de points sur l'axe X
NT = 100
dx = L/NX # distance entre les points
x = [] #Les différents points sur l'axe X
f = [] #les differentes valeurs de f pour chaque x
lapf = [] #laplacien
gradf = [] #gradient
divf = [] #divergence
n= random.randint(3, 200) #génère un int entre 3 et 200 compris
h= 1/(n-1)
normErrGrad = []
normErrLap = []

##FONCTIONS UTILES##
def norm2(list):
    norm2 = 0
    for i in range(len(list)):
        norm2 = norm2 + list[i]**2
    return np.sqrt(norm2)

##CREATION DES POINTS SUR L'AXE X##
x = np.arange(0, L, dx) #On crée tous les points entre 0 et 5 espacés régulièrement

##CREATION DE LA FONCTION##
def u(x):
    f=[]
    for i in range (0, NX):
        f = f + [(x[i])**4] #On construit la liste des différentes valeurs prises par la fonction, ici x**4
    return f

f=u(x)

##CALCUL DU GRADIENT##
for i in range(len(x)-2):
    gradf.append( (f[i]+f[i+2])/(2*dx) ) #marche que jusqu'à 97, car on fait x[i+2] et x a des coordonnées de 0 à 99
for i in range(len(x)-2, len(x)): #il faut créer les deux dernières coordonnées du vecteur grad
    gradf.append( 2*gradf[i-1]-gradf[i-2] ) #rajoute une coordonnées à la liste
    
def ux(x):
    fx=[]
    for i in range (0, NX):
        fx = fx + [4*(x[i])**3] #On construit la liste des différentes valeurs prises par la fonction, ici x**4
    return fx

fx=ux(x)

##CALCUL DE LA DIVERGENCE##
for i in range(0, NX-2):
    divf.append( gradf[i] )
for i in range(len(x)-2, len(x)):
    divf.append( 2*divf[i-1]-divf[i-2] )

##CALCUL DU LAPLACIEN##
for i in range (0, NX-2):
    lapf.append( (gradf[i]+gradf[i+2])/(2*dx))
for i in range(len(x)-2, len(x)):
    lapf.append( (2*lapf[i-1]-lapf[i-2]) )

def uxx(x):
    fxx=[]
    for i in range (0, NX):
        fxx = fxx + [12*(x[i])**2] #On construit la liste des différentes valeurs prises par la fonction, ici x**4
    return fxx

fxx=uxx(x)

##CALCUL DES ERREURS##
errGrad= abs(np.array(gradf)-np.array(fx)) #liste de liste des erreurs en chaque point du maillage pour la première dérivée
errLap= abs(np.array(lapf)-np.array(fxx)) #liste de liste des erreurs pour la deuxième dérivée
normErrGrad = h * norm2(errGrad)
normErrLap = h*norm2(errLap)


plt.figure(1)
plt.plot(x,f,label='courbe de f',color='black')
plt.plot(x,gradf,label='courbe de grad(f)',color='blue')
plt.plot(x,lapf,label='courbe de lapl(f)',color='green')
plt.legend((('f(x)'), ('grad(f)'), ('div(f)')))
plt.title('Projet 1')
plt.show()

plt.figure(2)
plt.plot(x,f,label='courbe de f',color='black')
plt.plot(x,gradf,label='courbe de grad(f)',color='blue')
plt.plot(x,lapf,label='courbe de lapl(f)',color='green')
plt.legend((('f(x)'), ('grad(f)'), ('div(f)')))
plt.title('Projet 1')
plt.show()









