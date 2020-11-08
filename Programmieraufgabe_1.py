#!/usr/bin/env python
# coding: utf-8

# # Programmieraufgabe 1: Gram-Schmidt Orthonormalisierungsverfahren

# In[1]:


import numpy as np
from scipy.linalg import norm, hilbert


# ## a)

# In[2]:
def skalarprodukt(V,W):
    n=len(V)
    x=0
    i=0
    while i<n:
        x=x+V[i]*W[i]
        i=i+1
    return x

def multskalar(a,W):
    n=len(W)
    B=[]
    i=0
    while i<n:
        B.append(a*W[i])
        i=i+1
    return B

def subtrahiere(V,W):
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
    print(n,m)
    Q=A #wird am ende ausgegeben
    spalte=1
    while spalte<m:
        print('spalte ',spalte+1)
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
            neuespalte=subtrahiere(neuespalte, C)
            i=i+1        
        t=0
        while t<n:
            Q[t][spalte]=neuespalte[t]
            t=t+1
        spalte=spalte+1
        
    return Q


# ## b)

# In[3]:


H = hilbert(30)

for n in [5,10,15,20]:
    ## Ihr Code hier


# ## c)

# In[4]:


def mod_gram_schmidt(V):
    Q = np.copy(V)
    ## Ihr Code hier
    
    return Q


# ## d)

# In[5]:


for n in [5,10,15,20]:
    ## Ihr Code hier


# In[ ]:




