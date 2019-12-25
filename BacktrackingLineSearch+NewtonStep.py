#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math


# In[2]:


def f_seq(X):
    return X.T@C@X

def f_hole(X):
    return 1-np.exp(-(X.T@C@X))

def df_seq(X):
    return 2*X.T@C

def df_hole(X):
    return np.exp(-(X.T@C@X))*2*C@X.T@C

def C_fn(entry, size):
    C_tmp=[math.pow(entry,(i-1)/(n-1)) for i in range(1,size+1)]
    return np.diag(C_tmp)

def GD_with_backtracking(f,df,X,alpha):
    count=0
    file= open("path.dat","w+")
    x_str=""
    for j in range(X.shape[0]):
        x_str+=str(X[j])+" "
    x_str+=str(f(X))+"\n"
    file.write(x_str)
    while(count<10):
        beta=-(df(X)/np.linalg.norm(df(X)))
        # beta=-df(X)       # Modification to see that the relation is not based on the mag. of gradient
        print("beta",beta)
        print("condition",(f(X+alpha*beta)> (f(X)+wolfe_para*df(X).T@(alpha*beta))).any())
        while (f(X+alpha*beta)> (f(X)+wolfe_para*df(X).T@(alpha*beta))).any():
            alpha=alpha_dec*alpha
            print("alpha",alpha)
        X=X+alpha*beta
        # Write in the file
        x_str=""
        for j in range(X.shape[0]):
            x_str+=str(X[j])+" "
        x_str+=str(f(X))+"\n"
        file.write(x_str)
        print("X",X)
        alpha=min(alpha_inc*alpha,delta_max)
        if (abs(alpha*beta) < tolerance_theta).all(): 
            print(count)
            count+=1
    file.close()
    return X,alpha

def GD_with_backtracking_newton(f,df,X,alpha):
    count=0
    file= open("path_n.dat","w+")
    x_str=""
    for j in range(X.shape[0]):
        x_str+=str(X[j])+" "
    x_str+=str(f(X))+"\n"
    file.write(x_str)
    while(count<10):
        beta=-(np.linalg.inv(C))@(df(X)/np.linalg.norm(df(X)))
        print("beta",beta)
        print("condition",(f(X+alpha*beta)> (f(X)+wolfe_para*df(alpha*beta).T)).any())
        while (f(X+alpha*beta)> (f(X)+wolfe_para*df(alpha*beta).T)).any():
            alpha=alpha_dec*alpha
            print("alpha",alpha)
        X=X+alpha*beta
        # Write in the file
        x_str=""
        for j in range(X.shape[0]):
            x_str+=str(X[j])+" "
        x_str+=str(f(X))+"\n"
        file.write(x_str)
        print("X",X)
        alpha=min(alpha_inc*alpha,delta_max)
        if (abs(alpha*beta) < tolerance_theta).all(): 
            print(count)
            count+=1
    file.close()
    return X,alpha


# In[3]:


#Initialize magic values
n=2
x0=[1. for i in range(n)]
c=10
tolerance_theta=0.1
alpha_inc=1.2
alpha_dec=0.5
delta_max=float("inf")
wolfe_para=0.01
alpha=1

# Get the n value according to the size of x0
X=np.array(x0)
print(X.shape)

# Initialize  C
C=C_fn(entry=c, size=n)


# In[4]:


# Test the f_seq(X)
print(df_hole(X))
# df_hole(X)/np.linalg.norm(df_hole(X))


# In[5]:


# Gradient Descent with Backtracking for the f_seq function
GD_with_backtracking(f=f_seq,df=df_seq,X=X,alpha=alpha)


# In[6]:


# Gradient Descent with Backtracking for the f_hole function
GD_with_backtracking(f=f_hole,df=df_hole,X=X,alpha=alpha)


# In[7]:


GD_with_backtracking_newton(f=f_seq,df=df_seq,X=X,alpha=alpha)

# GD_with_backtracking_newton(f=f_hole,df=df_hole,X=X,alpha=alpha)


# %%
