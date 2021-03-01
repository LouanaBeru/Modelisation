#Modules ÃƒÂ  importer
from math import sin,sqrt,exp,pi,log10
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#Initialisation 
k = 0.1
s= 200
Lx = 1
Nx = 200 #le maillage spatial sert a la representation de la solution et aux calculs des ps par integration numerique
hx = Lx/(Nx-1)
x = np.linspace(0,Lx,Nx)
modmax = 60
mode=[]
dffm=[]
#introduire une boucle sur modmax et calculer l'erreur ||u(modmax)||
f=np.zeros(len(x))
u0=np.zeros(len(x))
u=np.zeros(len(x))
testx=np.zeros(modmax)
fbx=np.zeros((modmax,Nx))      #phi_n
fmp=np.zeros(modmax)           #projection sur phi_n de f
imp=np.zeros(modmax)           #projection sur phi_n de u0
fun0=np.zeros(modmax)
epsi=[]


def norm2(x,y):
    norm2 = 0
    for i in range(modmax):
        norm2 = norm2 + (x[i]+y[i])**2
    return np.sqrt(norm2)


#rhs
#Coeur du programme

for i in range(1,Nx):
    f[i]=30*exp(-s*((x[i]-Lx/4)**2))  #(100*(x[i]**2)*((Lx-x[i])**2))/(Lx**4)
    u0[i]=exp(-s*((x[i]-Lx/2)**2)) #+exp(-2*s*((x[i]-Lx/3)**2))+exp(-3*s*((x[i]-2*Lx/3)**2))
    u[i]=u0[i]


plt.figure(1)
plt.plot(x,u0,'g')

espf= log10(norm2(fun0,u0)*2.7)


for m in range(1,modmax):
    epsi+=[espf] 

#fb normalise
# 
for m in range(1,modmax):
    for i in range(1,Nx):
        fbx[m][i]=sin(pi*x[i]/Lx*m)               #phi_n
        testx[m]=testx[m]+fbx[m][i]*fbx[m][i]*hx  
    testx[m]=sqrt(testx[m])                       #norme L2 de phi_n
    for i in range(1,Nx):                         #normalisation
        fbx[m][i]=fbx[m][i]/testx[m]

        
#verifier l'orthonormalite des fbx ?
#<fbx[m],fbx[n]>=delta_mn


#projection f second membre et u0 condi init sur fbx

for m in range(1,modmax):
    for i in range(1,Nx):
        fmp[m]+=f[i]*fbx[m][i]*hx    # <f,\phi_n> = f_n 
        imp[m]+=u0[i]*fbx[m][i]*hx   # <u0,phi_n> = c_n

#somme serie

temps=0.0    #on doit retrouver la condition initiale
for i in range(1,Nx):
    u[i]=0
for m in range(1,modmax):
    al=(m**2)*(pi**2)/(Lx**2)*k
    coef=imp[m]*exp(-al*temps)
    for i in range(0,Nx-1):
         u[i]+=fbx[m][i]*coef
    dffm += [norm2(u0,u)]
    mode+=[m]     
 
plt.plot(x,u,'blue')
 
temps=0.02  #la solution a n'importe quel temps sans avoir a calculer les iter intermediaires
for i in range(1,Nx):
    u[i]=0
for m in range(1,modmax):
    al=(m**2)*(pi**2)/(Lx**2)*k
    coef=imp[m]*exp(-al*temps)
    coeff=fmp[m]*(1-exp(-al*temps))/al
    for i in range(0,Nx):
        u[i]+=fbx[m][i]*(coeff+coef)

plt.plot(x,u,'r')

plt.figure(2)


for m in range(modmax-1):
    dffm[m]= log10(dffm[m])
plt.plot(mode,dffm,'r')
plt.plot(mode,epsi,'g')

plt.show()
print(dffm[23])
print(espf)

while dffm[23]>espf:
    print(Nx)
    print(log10(dffm[23]))
    print(espf)
    Nx -= 1 #le maillage spatial sert a la representation de la solution et aux calculs des ps par integration numerique
    x = np.linspace(0,Lx,Nx)
    hx = Lx/(Nx-1)
    mode=[]
    dffm=[]
    #introduire une boucle sur modmax et calculer l'erreur ||u(modmax)||
    f=np.zeros(len(x))
    u0=np.zeros(len(x))
    u=np.zeros(len(x))
    testx=np.zeros(modmax)
    fbx=np.zeros((modmax,Nx))      #phi_n
    fmp=np.zeros(modmax)           #projection sur phi_n de f
    imp=np.zeros(modmax)           #projection sur phi_n de u0
    fun0=np.zeros(modmax)



    #rhs
    #Coeur du programme

    for i in range(1,Nx):
        f[i]=30*exp(-s*((x[i]-Lx/4)**2))  #(100*(x[i]**2)*((Lx-x[i])**2))/(Lx**4)
        u0[i]=exp(-s*((x[i]-Lx/2)**2)) #+exp(-2*s*((x[i]-Lx/3)**2))+exp(-3*s*((x[i]-2*Lx/3)**2))
        u[i]=u0[i]




    #fb normalise
    # 
    for m in range(1,modmax):
        for i in range(1,Nx):
            fbx[m][i]=sin(pi*x[i]/Lx*m)               #phi_n
            testx[m]=testx[m]+fbx[m][i]*fbx[m][i]*hx  
        testx[m]=sqrt(testx[m])                       #norme L2 de phi_n
        for i in range(1,Nx):                         #normalisation
            fbx[m][i]=fbx[m][i]/testx[m]

        
    #verifier l'orthonormalite des fbx ?
    #<fbx[m],fbx[n]>=delta_mn


    #projection f second membre et u0 condi init sur fbx

    for m in range(1,modmax):
        for i in range(1,Nx):
            fmp[m]+=f[i]*fbx[m][i]*hx    # <f,\phi_n> = f_n 
            imp[m]+=u0[i]*fbx[m][i]*hx   # <u0,phi_n> = c_n

    #somme serie

    temps=0.0    #on doit retrouver la condition initiale
    for i in range(1,Nx):
        u[i]=0
    for m in range(1,modmax):
        al=(m**2)*(pi**2)/(Lx**2)*k
        coef=imp[m]*exp(-al*temps)
        for i in range(0,Nx-1):
             u[i]+=fbx[m][i]*coef
        dffm += [norm2(u0,u)]
        mode+=[m]     

 
    temps=0.02  #la solution a n'importe quel temps sans avoir a calculer les iter intermediaires
    for i in range(1,Nx):
        u[i]=0
    for m in range(1,modmax):
        al=(m**2)*(pi**2)/(Lx**2)*k
        coef=imp[m]*exp(-al*temps)
        coeff=fmp[m]*(1-exp(-al*temps))/al
        for i in range(0,Nx):
            u[i]+=fbx[m][i]*(coeff+coef)


print(Nx+1)


