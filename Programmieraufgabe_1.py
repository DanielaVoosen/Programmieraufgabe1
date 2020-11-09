#!/usr/bin/env python
# coding: utf-8

# # Programmieraufgabe 1: Gram-Schmidt Orthonormalisierungsverfahren

# In[1]:


import numpy as np
from scipy.linalg import norm, hilbert


# ## a)

# In[2]:
def skalarprodukt(V,W):#Normales Skalarprodukt
    n=len(V)
    x=0
    i=0
    while i<n:
        x=x+V[i]*W[i]
        i=i+1
    return x

def multskalar(a,W):#Multiplikation: Skalar, Vektor
    n=len(W)
    B=[]
    i=0
    while i<n:
        B.append(a*W[i])
        i=i+1
    return B

def subtrahiere(V,W):#Vektorsubtraktion
    n=len(V)
    B=[]
    i=0
    while i<n:
        B.append(V[i]-W[i])
        i=i+1
    return B

def gram_schmidt(A):
    n=len(A)
    m=len(A[0])
    Q=A #wird am ende ausgegeben
    spalte=1# da Q=A kann in Spalte 1, bzw 2, angefangen werden
    while spalte<m:# die restlichenw erden durchgegangen
        altespalte=[]
        z=0
        V=[]
        while z<n:
                V.append(Q[z][spalte-1])#V entpricht v_i
                altespalte.append(A[z][spalte])#erstes Element aus Iteration
                z=z+1
        i=0
        neuespalte=altespalte
        while i<spalte:
            a=skalarprodukt(V,altespalte)
            b=skalarprodukt(V,V)
            C=multskalar(a/b,V)
            neuespalte=subtrahiere(neuespalte, C)#entspricht der Summe
            i=i+1        
        t=0
        while t<n:
            Q[t][spalte]=neuespalte[t]#Eintr?ge in Q werden ge?ndert
            t=t+1
        spalte=spalte+1
    return Q
        
        #Testmatrix:
A=[[3,2],[1,2],[2,2]]
B=gram_schmidt(A)
print(B)
#Ergebnnis nach Wikipedia: [[3,-4/7],[1,8/7],[2,2/7]]

    return Q





# ## b)

H = hilbert(30)

for n in [5,10,15,20]:
   i=29 #da zählen bei 0 anfängt
    H_2=H
    m=len(H_2)
    while i>=n:#Es müssen 30-n Spalten entfernt werden
        h=0
        while h<m:
            del(H_2[h][i])#in jeder Zeile wird das ite Element entfernt
            h+=1
        i-=1
    Q=gram_schmidt(H_2)
    print('Für n= ',n,' gibt das Gram-Schmidt-Verfahren die Matrix ',Q,' zurück')
    
    #HIER FEHLT DER ZWEITE TEIL


# ## c)

# In[4]:


def mod_gram_schmidt(V):
    Q = np.copy(V)
    ## Ihr Code hier
    
    return Q


# ## d)

# In[5]:
H = hilbert(30)

for n in [5,10,15,20]:
   i=29 #da zählen bei 0 anfängt
   H_2=H
   m=len(H_2)
   while i>=n:#Es müssen 30-n Spalten entfernt werden
        h=0
        while h<m:
            del(H_2[h][i])#in jeder Zeile wird das ite Element entfernt
            h+=1
        i-=1
   Q=gram_schmidt(H_2)
   print('Für n= ',n,' gibt das Gram-Schmidt-Verfahren die Matrix ',Q,' zurück')


# In[ ]:




