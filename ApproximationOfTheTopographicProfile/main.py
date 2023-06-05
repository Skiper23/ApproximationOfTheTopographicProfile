import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bisect
WIELOMIAN = 'w'
POCHODNA = 'p'
DRUGAPOCHODNA = 'd'
DRUGAPOCHODNAZERO= 'z'
def iloczynMianownik(indeks, wsp):
    sum=1
    for j in range(0,len(wsp)):
        if indeks==j:
            continue
        sum=sum*(wsp[indeks]-wsp[j])
    return sum

def iloczynLicznik(x,wsp,skip):
    sum=1
    for i in range(0,len(wsp)):
        if i==skip:
            continue
        sum=(x-wsp[i])*sum
    return sum
def wartoscInterpolacji(bazaLagrange, wspy):
    sum=0
    for i in range(0, len(wspy)):
        sum=sum+bazaLagrange[i]*wspy[i]
    return sum
def interpolacjaLagrange(xdata,x):
    for i in xdata:
        for j in range(0, len(x)):
            l = iloczynLicznik(i, x, j)
            m = iloczynMianownik(j, x)
            bazaLagrange.append(l / m)
        ywyliczone.append(wartoscInterpolacji(bazaLagrange, y))
        bazaLagrange.clear()

    plt.subplot()
    plt.plot(xdata, ywyliczone, color="blue")
    plt.plot(xdata, ydata, color="red")
    plt.show()

def tworzWspolczynniki(xn,x,indeks,typ,size):
    wektor=np.zeros(size)

    if typ == WIELOMIAN:
        wektor[4 * indeks] = 1
        if xn == x:
            return wektor
        wektor[4 * indeks + 1]=xn-x
        wektor[4 * indeks + 2] = pow(xn - x,2)
        wektor[4 * indeks + 3] = pow(xn - x,3)
    elif typ == POCHODNA:
        wektor[4 * indeks + 1]=1
        wektor[4 * indeks + 2] = (xn - x)*2
        wektor[4 * indeks + 3] = pow(xn - x,2)*3
        wektor[4 * (indeks+1) + 1] = -1
    elif typ == DRUGAPOCHODNA:
        wektor[4 * indeks + 2] = 2
        wektor[4 * indeks + 3] = (xn - x)*6
        wektor[4 * (indeks+1) + 2] = -2
    elif typ == DRUGAPOCHODNAZERO:
        wektor[4 * indeks + 2] = 2
        wektor[4 * indeks + 3] = (xn - x)*6
    return wektor
def tworzMacierz(xdata,ydata):
    size=4*(len(xdata)-1)
    macierz=np.zeros((size,size))
    wektor=np.zeros((size,1))
    indeks=0

    macierz[indeks, :] = tworzWspolczynniki(xdata[0], xdata[0], 0, WIELOMIAN,size)
    wektor[indeks] = ydata[0]
    indeks+=1
    macierz[indeks, :] = tworzWspolczynniki(xdata[1], xdata[0], 0, WIELOMIAN,size)
    wektor[indeks] = ydata[1]
    indeks += 1
    for i in range(1,(len(xdata)-1)):
        macierz[indeks,:]=tworzWspolczynniki(xdata[i], xdata[i], i, WIELOMIAN,size)
        wektor[indeks]=ydata[i]
        indeks+=1
        macierz[indeks,:]=tworzWspolczynniki(xdata[i+1], xdata[i], i, WIELOMIAN,size)
        wektor[indeks]=ydata[i+1]
        indeks+=1
        macierz[indeks,:]=tworzWspolczynniki(xdata[i+1], xdata[i], i-1, POCHODNA,size)
        wektor[indeks]=0
        indeks+=1
        macierz[indeks,:]=tworzWspolczynniki(xdata[i+1], xdata[i], i-1, DRUGAPOCHODNA,size)
        wektor[indeks]=0
        indeks+=1

    macierz[indeks, :] = tworzWspolczynniki(xdata[0], xdata[0], 0, DRUGAPOCHODNAZERO,size)
    wektor[indeks] = 0
    indeks += 1
    macierz[indeks, :] = tworzWspolczynniki(xdata[len(xdata)-1], xdata[len(xdata)-2], len(xdata)-2, DRUGAPOCHODNAZERO,size)
    wektor[indeks] = 0
    indeks += 1
    return np.linalg.solve(macierz,wektor)

def sklejWartosci(wektorRozwiazan, xdane, x):
    indeks=bisect.bisect(xdane,x)
    sum=0
    for i in range(4*(indeks-1), 4*(indeks-1)+4):
        pom=wektorRozwiazan[i]*pow(x-xdane[indeks-1],i-4*(indeks-1))
        sum=sum+pom
    return sum


dane = pd.read_csv("2018_paths//Obiadek.csv")
xdata = dane['x'].tolist()
ydata = dane['y'].tolist()
x=xdata[::10]
y=ydata[::10]
xtest=xdata[::10]
ytest=ydata[::10]
bazaLagrange=[]
ywyliczone=[]

#interpolacjaLagrange(xdata,x)
wektorRozwiazan=tworzMacierz(x,y)
print(len(x))
print(len(wektorRozwiazan))
for i in xdata:
    if(i<x[len(x)-1]):
        ywyliczone.append(sklejWartosci(wektorRozwiazan, x, i))
    else:
        ywyliczone.append(0)

plt.subplot()
plt.plot(xdata, ywyliczone, color="blue")
plt.plot(xdata, ydata, color="red")
plt.show()
print(ywyliczone)
print(ydata)
print(x)
print(y)