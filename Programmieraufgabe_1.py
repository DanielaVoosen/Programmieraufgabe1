#!/usr/bin/env python
# coding: utf-8

# # Programmieraufgabe 1: Gram-Schmidt Orthonormalisierungsverfahren

# In[1]:
import math
import numpy as np
from scipy.linalg import norm, hilbert


# ## a)

# In[2]:
def skalarprodukt(V,W):#Normales Skalarprodukt
    n = len(V)
    x = 0
    i = 0
    while i < n:
        x = x + V[i] * W[i]
        i = i + 1
    return x

def multskalar(a,W):#Multiplikation: Skalar, Vektor
    n = len(W)
    B = []
    i = 0
    while i < n:
        B.append(a * W[i])
        i = i + 1
    return B

def subtrahiere(V,W):#Vektorsubtraktion
    n = len(V)
    B = []
    i = 0
    while i < n:
        B.append(V[i] - W[i])
        i = i + 1
    return B

def gram_schmidt(A):
    n = len(A)
    m = len(A[0])
    Q = A #wird am Ende ausgegeben
    spalte = 1# da Q=A kann in Spalte 1, bzw 2, angefangen werden
    while spalte < m:# die restlichen werden durchgegangen
        altespalte = []
        z = 0
        V = []
        while z < n:
                V.append(Q[z][spalte - 1])#V entpricht v_i
                altespalte.append(A[z][spalte])#erstes Element aus Iteration
                z = z + 1
        i = 0
        neuespalte = altespalte
        while i < spalte:
            a = skalarprodukt(V,altespalte)
            b = skalarprodukt(V,V)
            C = multskalar(a / b,V)
            neuespalte = subtrahiere(neuespalte, C)#entspricht der Summe
            i = i + 1        
        t = 0
        while t < n:
            Q[t][spalte] = neuespalte[t]#Eintraege in Q werden geaendert
            t = t + 1
        spalte = spalte + 1
    return Q
    
#Testmatrix:
A = [[3,2],[1,2],[2,2]]
B = gram_schmidt(A)
print(B)
#Ergebnnis nach Wikipedia: [[3,-4/7],[1,8/7],[2,2/7]]



# ## b)
H = hilbert(30)

for n in [5,10,15,20]: #wir bauen uns hier eine kleinere Hilbert-Matrix (also aus einer 30x30 Matrix
                       #wird eine 30x5, 30x10,,,)

    Q_GS = gram_schmidt(H[0:n,0:n]) #Das ist die Macht von Numpy :D Vorteil: Zeilenumbrueche werden richtig gesetzt
    print('Für n =',n,' gibt das Gram-Schmidt-Verfahren folgende Matrix zurueck: \n\n\n', Q_GS, "\n\n") 

    #m=len(H) # m: Anzahl der Zeilen, n: Anzahl der Spalten
    #j=0
    #HB=[]
    #while j<m:
    #    V=[]
    #    i=0
    #    while i<n:
    #        V.append(H[j][i])
    #        i+=1
    #    HB.append(V)
    #    j+=1
    #Q=gram_schmidt(HB)
    #print('Für n= ',n,' gibt das Gram-Schmidt-Verfahren folgende Matrix
    #zurueck: \n', Q)

    I = np.identity(n) #die Einheitsmatrix mit Laenge/Breite n wird gebaut
    Q_GS_T = np.transpose(Q_GS) #wir transponieren die Matrix Q
   
    pruefen_GS = (1 / n) * norm(I - np.matmul(Q_GS_T,Q_GS)) #wir pruefen, ob die Spalten von Q ein Orthonormalsystem bilden (Formel vom
                                                   #Blatt)
    print('Pruefung durch gram_schmidt: n^(-1)||I-Q^(T)Q|| = ',pruefen_GS,"\n\n") #was faellt uns auf? Die Werte werden immer kleiner, desto groesser die Matrix?
    #Desto hoeher das n, desto mehr ist das System ein Orthogonalsystem.
    #Desto hoeher das n, desto kleiner werden die Werte von links oben nach rechts unten in der Hilbertmatrix


# ## c)

# In[4]:
def mod_gram_schmidt(A):
    Q = A #wie in der a
    n = len(A)
    i = 1
    while i < n: #wir gehen wieder alle Spalten durch
        U = []
        a = 0
        while a < n: #startvektor der Iteration
            U.append(Q[a][i])
            a+=1
        j = 0
        while j < i: #hier wird U geändert, wie in der Formel
           z = skalarprodukt(U,Q[i - 1])
           U = subtrahiere(U,multskalar(z,Q[i - 1]))
           j+=1
        norm = math.sqrt(skalarprodukt(U,U))#letztendlich noch deren Norm
        d = 0
        p = 1 / norm
        X = multskalar(p,U)#das ist Definition der neuen Spalte
        while d < n:# Uebertragen der Spalte auf die Ausgabematrix
            Q[d][i] = X[d]
            d+=1
        i+=1
    return Q
    

# ## d)

# In[5]:
for n in [5,10,15,20]:
    Q_MGS = mod_gram_schmidt(H[0:n,0:n]) #Das ist die Macht von Numpy :D Vorteil: Zeilenumbrueche werden richtig gesetzt
    print('Für n =',n,' gibt das Gram-Schmidt-Verfahren folgende Matrix zurueck: \n\n\n', Q_MGS, "\n\n") 
    I = np.identity(n) #die Einheitsmatrix mit Laenge/Breite n wird gebaut
    Q_MGS_T = np.transpose(Q_MGS) #wir transponieren die Matrix Q
   
    pruefen_MGS = (1 / n) * norm(I - np.matmul(Q_MGS_T,Q_MGS)) #wir pruefen, ob die Spalten von Q ein Orthonormalsystem bilden (Formel vom
                                                   #Blatt)
    print('Pruefung durch mod_gram_schmidt: n^(-1)||I-Q^(T)Q|| = ',pruefen_MGS,"\n\n") 

#for n in [5,10,15,20]:
#    m = len(H)
#    j = 0
#    HB = []
#    while j < m:
#        V = []
#        i = 0
#        while i < n:
#            V.append(H[j][i])
#            i+=1
#        HB.append(V)
#        j+=1
#    Q = mod_gram_schmidt(HB)
#    print('Für n= ',n,' gibt das modifizierte Gram-Schmidt-Verfahren die Matrix ',Q,' zurück')


# In[ ]:




