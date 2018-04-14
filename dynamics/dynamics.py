
# coding utf-8

# In[2]:


import pandas as pd
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


def IDlimb(Ax,Ay,MA,theta,alpha,m,I,l,lc,ax,ay):
    Kx=np.zeros(len(theta))
    Ky=np.zeros(len(theta))
    Mk=np.zeros(len(theta))
    for i in range(1,len(theta)):
        Kx[i]=Ax[i]+m*ax[i];
        Ky[i]=Ay[i]+m*ay[i]+m*9.81;
        Mk[i]=MA[i]+Ax[i]*l*math.sin(theta[i])-Ay[i]*l*math.cos(theta[i])-m*9.81*lc*math.cos(theta[i])+(I+m*lc*lc)*alpha[i];
    return Kx,Ky,Mk


# In[3]:


def IDfoot(F1x,F1y,theta,aplha,m,I,l,lc,ax,ay):
    F2x=np.zeros(len(theta))
    F2y=np.zeros(len(theta))
    M2=np.zeros(len(theta))
    for i in range(1,len(theta)):
        F2x[i]=m*ax[i]-F1x[1]
        F2y[i]=m*9.81+m*ay[1]-F1y[1]
        M2[i]=(I+m*lc*lc)*aplha[i]+F1y[i]*math.cos(theta[i])*l-F1x[i]*math.sin(theta[i])*l+m*9.81*math.cos(theta[i])*(l-lc)
    return F2x,F2y,M2


# In[6]:


def jointangle(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,n):
    for i in range(0,n):
        m1[i]=(pt2y[i]-pt1y[i])/(pt2x[i]-pt1x[i])
        m2[i]=(pt3y[i]-pt1y[i])/(pt3x[i]-pt1x[i])
        angle[i]=np.arctan((m2[i]-m1[i])/(1+m1[i]*m2[i]))
    return angle


# In[8]:


def angularaccel(w1,n,sampling_time):         #j1 is joint angle dataframe #n is number of datapoints 
    for i in range(1,n-2):
        ang_vel[i-1] = (w1[i+1,0]-w1[i-1,0])/(2* sampling_time)
    return ang_vel


# In[9]:


def derivative(j1,time): #j1 is joint angle dataframe #n is number of datapoints 
    deriv=np.zeros(len(j1))
    for i in range(1,len(j1)-1):
        deriv[i] = (j1[i+1]-j1[i-1])/(time[i+1]-time[i-1])
    return deriv


# In[13]:


def acceleration(l,lc,alpha,omega,theta):
    r=l-lc
    ax=np.zeros(len(alpha))
    ay=np.zeros(len(alpha))
    for i in range(1,len(alpha)):
        ax[i]=-(r*alpha[i]*math.sin(theta[i])+omega[i]*omega[i]*r*math.cos(theta[i]))
        ay[i]=r*alpha[i]*math.cos(theta[i])-omega[i]*omega[i]*r*math.sin(theta[i])
    return ax,ay


# In[4]:


def absangle(pt1x,pt1y,pt2x,pt2y,n):
    angle = np.zeros(len(pt1x))
    for i in range(0,len(pt1x)):
        angle[i]=np.arctan2((pt2y[i]-pt1y[i]),(pt2x[i]-pt1x[i]))
    return angle


# In[12]:


def UpBod(F1x,F1y,M1,I,l,theta,mB,n):
    for i in range(1,n):    
        ax[i]=-F1x[i]/mB
        ay[i]=-(mB*9.81+F1y[i])/mB
        alpha[i]=-M1[i]/(I+mB*l*l)
    return ax,ay,alpha

def AnthroData(weight,height):
    data=np.zeros(15)
    misc= (7.8*9.6*9.6*49.5)+(46.84*31.6*31.6*50.3)+2*(2.7*16.4*16.4*32.3)+2*(2.7*13.7*13.7*30.3)+2*(0.6*8.2*8.2*29.7) #sum of Icm of all the upper body parts
    data[0]= weight*67.2/100 #weight of UB
    data[1]= height*14.2/100 #height of com of UB from hip 
    data[2]= data[0]*data[1]*data[1]+misc #moment for UB
    data[3]= weight*9.9/100#weight of thigh
    data[4]= height*25.4/100#height of thigh
    data[5]= data[4]*43.3/100#height of com thigh
    data[6]= data[3]*data[4]*data[4]*32.3/100#moment of thigh
    data[7]= weight*4.6/100#weight of shank
    data[8]= height*23.3/100#height of shank
    data[9]= data[8]*43/100 #height of com shank
    data[10]= data[7]*data[8]*data[8]*30.2/100 #moment of shank
    data[11]= weight*1.4/100#weight of foot
    data[12]= height*11.7/100 #height of foot
    data[13]= data[12]*50/100 #height of com foot
    data[14]= data[11]*data[12]*data[12]*47.5/100 #moment of foot
    return data