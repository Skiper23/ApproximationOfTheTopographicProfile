import pandas as pd
import copy
def horner(x, j):
    sum=1
    for i in range(0,len(x)):
        if i==j:
            continue
        sum=sum*(x[j]-x[i])
    return sum
def sumaWspolczynnikow(x, pom):
    for i in range(0,len(pom)-1):
        x[i+1]=x[i+1]+pom[i]
    return x

def bazaLagrange(x,k):
    wspolczynniki=[1]
    pom=wspolczynniki.copy()
    for i in range(0,len(x)):
        if i==k:
            continue
        pom=wspolczynniki.copy()
        for j in range(len(wspolczynniki)):
            wspolczynniki[j]=wspolczynniki[j]*(-1)*x[i]
        wspolczynniki.append(1)
        wspolczynniki=sumaWspolczynnikow(wspolczynniki,pom)

    return wspolczynniki

dane = pd.read_csv("2018_paths//100.csv")
x = dane['x'].tolist()
y = dane['y'].tolist()
x=x[:82]
y=y[:82]
wsplist=[]
wspmianownik=[]
for i in range(len(x)):
    wspmianownik.append(horner(x,i))
    wsplist.append(bazaLagrange(x,i))

for i in range(len(wsplist)):
    for j in range(len(wsplist[i])):
        wsplist[i][j]=wsplist[i][j]/wspmianownik[i]*y[i]

wyniki=[]
for i in range(len(wsplist)):
    sum=0
    for j in range(len(wsplist)):
        sum=sum+wsplist[j][i]
    wyniki.append(sum)
print(wyniki)