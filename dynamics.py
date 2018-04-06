
# coding: utf-8

# In[5]:


import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# In[1]:


def filterdata(x,y,n):
    from scipy.signal import filtfilt
    # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    yy = filtfilt(b,a,y)
    plt.plot(x, yy, linewidth=2, linestyle="-", c="b")
    plt.plot(x,y,linewidth=2, linestyle="-", c="r")
    return yy;


# In[2]:


def IDlimb(Ax,Ay,MA,theta_s,alpha_s,ms,Is,ls,lcs,ax,ay,n):
    for i in range(1,n):
        Kx[i]=Ax(i)+ms*ax(i);
        Ky[i]=Ay(i)+ms*ay(i)+ms*9.81;
        Mk[i]=MA(i)+Ax(i)*ls*sin(theta_s(i))-Ay(i)*ls*Cos(theta_s(i))-ms*9.81*lcs*Cos(theta_s(i))+(Is+ms*lcs^2)*alpha_s(i);
    return Kx,Ky,Mk


# In[3]:


def IDfoot(F1x,F1y,theta,aplha,ax,ay,m,l,lc,I,n):
    for i in range(1,n):
        F2x[i]=m*ax[i]-F1x[1]
        F2y[i]=m*9.81+m*ay[1]-F1y[1]
        M2[i]=(I+m*lc*lc)*aplha[i]+F1y[i]*cos(theta[i])*l-F1x[i]*sin(theta[i])*l+m*9.81*cos(theta[i])*(l-lc)
    return F2x,F2y,M2


# In[6]:


def jointangle(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,n):
    for i in range(0,n):
        m1[i]=(pt2y[i]-pt1y[i])/(pt2x[i]-pt1x[i])
        m2[i]=(pt3y[i]-pt1y[i])/(pt3x[i]-pt1x[i])
        angle[i]=atan((m2[i]-m1[i])/(1+m1[i]*m2[i]))
    return angle


# In[8]:


def angularaccel(w1,n,sampling_time):         #j1 is joint angle dataframe #n is number of datapoints 
    for i in range(1,n-2):
        ang_vel[i-1] = (w1[i+1,0]-w1[i-1,0])/(2* sampling_time)
    return ang_vel


# In[9]:


def angularvel(j1,n,sampling_time):         #j1 is joint angle dataframe #n is number of datapoints 
    for i in range(1,n-2):
        ang_vel[i-1] = (j1[i+1,0]-j1[i-1,0])/(2* sampling_time)
    return ang_vel


# In[13]:


def acceleration(l,lc,aplha,omega,theta,n):
    r=l-lc
    for i in range(1,n):
        ax[i]=-(r*alpha[i]*sin(theta[i])+omega[i]*omega[i]*r*cos(theta[i]))
        ay[i]=r*alpha[i]*cos(theta[i])-omega[i]*omega[i]*r*sin(theta[i])
    return ax,ay


# In[11]:


def absangle(pt1x,pt1y,pt2x,pt2y,n):
    for i in range(0,n):
        angle[i]=atan(pt2y[i]-pt1y[i])/(pt2x[i]-pt1x[i])
    return angle


# In[12]:


def UpBod(F1x,F1y,M1,I,l,theta,mB,n):
    for i in range(1,n):    
        ax[i]=-F1x[i]/mB
        ay[i]=-(mB*9.81+F1y[i])/mB
        alpha[i]=-M1[i]/(I+mB*l*l)
    return ax,ay,alpha

