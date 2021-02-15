import matplotlib.pyplot as plt 
import numpy as np
import math

#Fonction : Gradient en dimension 1
#Variables : func : la fonction dont on veut le gradient, [xmin,xmax] : l'intervalle, NX : nbr de pts du maillage, epsdf : précision
#Renvoie : La liste du gradient de func évalué en tous les points du maillage

def Gradient(func,xmin, xmax, NX, epsdf):
    
    x=np.linspace(xmin, xmax, NX)
    dfdx = np.zeros(NX)

    for i in range(NX-1):
        dfdx[i]=(func(x[i]+epsdf)-func(x[i]-epsdf))/(2*epsdf)
    
    dfdx[NX-1] = 2*dfdx[NX-2] - dfdx[NX-3]

    return(dfdx)

#Fonction : Minimum en dimension 1
#Variables : Gradient(func), [xmin,xmax] : l'intervalle, NX : nbr de pts du maillage, epsdf : précision du minimum
#Renvoie : le minimum de func

def Minimum(func,xmin, xmax, NX, epsdf):

    dfdx = Gradient(func, xmin, xmax, NX, epsdf)
    x = np.linspace(xmin, xmax, NX)
    listMin = []

    localMinIndex = x[0]
    localMin = func(x[0])
    for i in range(NX-1):
        if (dfdx[i]*dfdx[i+1] <= 0 ):
            val = func(x[i])
            if val < localMin:
                localMin = val
                localMinIndex = x[i]
    return(localMinIndex)

##TEST avec la fonction f(x) = x**2
  
# def func(x):
#     f = math.sin((x+2)**2) - math.cos(x+2) 
#     return f

# x= np.linspace(-4, 4, 400)
# f = np.zeros(len(x))
# for i in range(len(x)):
#     f[i]=func(x[i])
# dfdx = Gradient(func, -4, 4, 400, 0.0001)

# plt.figure()
# plt.plot(x, f, color = 'blue')
# plt.plot(x, dfdx, color = 'yellow')
# m = Minimum(func, -4, 4, 400, 0.0001)
# plt.scatter(m, func(m), marker = 'P', color = 'green')
# plt.show()

## Que faire avec ce dfdx0, qu'est-ce qu'il représente, pourquoi a-t-on besoin de lui pour calculer le minimum de la fonction?
## Régler ça ( la fonction doit trouver le minimum) et ensuite c'est bon, il faudra l'appeler dans le programme projet3 en faisant import nomDuFichier.fonction()
## Je pourrais l'utiliser pour calculer le gradient de T et ensuite calculer le minimum de J et le tour est gagné !
## Pour calculer le minimum, j'avance, jusqu'à ce que le gradient de la fonction s'annule -> tester si la fonction est croissante
## Tester si c'est un minimum ou un maximum, essayer de renvoyer tous les minimums et à la fin de faire min, pour avoir le min global.
## Je dois recalculer le gradient pour calculer le minimum, et je dois faire un pas qui bouge ( plus grand si ça descend beaucoup, pus faible si on s'approche du min)
