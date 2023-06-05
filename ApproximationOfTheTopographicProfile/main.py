import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bisect
WIELOMIAN = 'w'
POCHODNA = 'p'
DRUGAPOCHODNA = 'd'
DRUGAPOCHODNAZERO= 'z'
ILOSCWEZLOW=30
FILE="2018_paths//WielkiKanionKolorado.csv"
def wyswietlWykres(xdata,ydata,ywyliczone,wezelx,wezely):
    plt.subplot()
    plt.plot(xdata, ywyliczone, color="red", label="")
    plt.ylim(500, 2900)
    plt.plot(xdata, ydata, color="blue")
    plt.scatter(wezelx, wezely, color="green")
    plt.title("Wielki Kanion Kolorado - ilosc wezlow 30")
    plt.xlabel("Dystans (m)")
    plt.ylabel("Wysokosc (m)")
    plt.legend(['F(x)','f(x)','f(x0)'])

    plt.show()

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
    bazaLagrange = []
    ywyliczone = []
    for i in xdata:
        for j in range(0, len(x)):
            l = iloczynLicznik(i, x, j)
            m = iloczynMianownik(j, x)
            bazaLagrange.append(l / m)
        ywyliczone.append(wartoscInterpolacji(bazaLagrange, wezely))
        bazaLagrange.clear()

    wyswietlWykres(xdata,ydata,ywyliczone,wezelx,wezely)

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
    if indeks == len(xdane):
        indeks-=1
    sum=0
    for i in range(4*(indeks-1), 4*(indeks-1)+4):
        pom=wektorRozwiazan[i]*pow(x-xdane[indeks-1],i-4*(indeks-1))
        sum=sum+pom
    return sum

def interpolacjafunkcjamisklejanymi(wezelx, wezely):
    wektorRozwiazan = tworzMacierz(wezelx, wezely)
    ywyliczone=[]
    for i in xdata:
        wartosc = sklejWartosci(wektorRozwiazan, wezelx, i)
        ywyliczone.append(float(wartosc))

    wyswietlWykres(xdata,ydata,ywyliczone,wezelx,wezely)


dane = pd.read_csv(FILE)
xdata = dane['x'].tolist()
ydata = dane['y'].tolist()
krok=round(len(xdata)/ILOSCWEZLOW)
wezelx=xdata[::krok]
wezely=ydata[::krok]


interpolacjaLagrange(xdata,wezelx)
interpolacjafunkcjamisklejanymi(wezelx, wezely)
