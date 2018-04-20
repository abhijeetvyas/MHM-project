import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

class BodyPart:
    time=[]
    def __init__(self, length,mass, CoM, Inertia,distal,proximal):
        self.sample=len(distal[0])
        self.length=length
        self.CoM=CoM
        self.mass=mass
        self.Inertia=Inertia
        self.distal=distal
        self.proximal=proximal
        self.ProximForce=[np.zeros(self.sample),np.zeros(self.sample)]
        self.ProximMoment=np.zeros(self.sample)
        self.angle = np.zeros(self.sample)
        self.omega=np.zeros(self.sample)
        self.alpha=np.zeros(self.sample)
        self.velProxim=[np.zeros(self.sample),np.zeros(self.sample)]
        self.accelProxim=[np.zeros(self.sample),np.zeros(self.sample)]
        self.accelCoM=[np.zeros(self.sample),np.zeros(self.sample)]
    def absangle(self):
        unfil=np.zeros(self.sample)
        for i in range(0,self.sample):
            unfil[i]=np.arctan2((self.proximal[1][i]-self.distal[1][i]),(self.proximal[0][i]-self.distal[0][i]))
        self.angle=filterdata(unfil,5)
        return self.angle

    def omega(self): #j1 is joint angle dataframe #n is number of datapoints
        unfil=np.zeros(self.sample)
        BodyPart.absangle(self)
        for i in range(1,self.sample-1):
            unfil[i] = (self.angle[i+1]-self.angle[i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
            self.omega=filterdata(unfil,5)
        return self.omega

    def alpha(self): #j1 is joint angle dataframe #n is number of datapoints 
        BodyPart.omega(self)
        unfil=np.zeros(self.sample)
        for i in range(1,self.sample-1):
            unfil[i] = (self.omega[i+1]-self.omega[i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
            self.alpha=filterdata(unfil,5)
        return self.alpha

    def velProxim(self): #j1 is joint angle dataframe #n is number of datapoints 
        unfil=[np.zeros(self.sample),np.zeros(self.sample)]
        for i in range(1,self.sample-1):
            unfil[0][i] = (self.proximal[0][i+1]-self.proximal[0][i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
            unfil[1][i] = (self.proximal[1][i+1]-self.proximal[1][i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
        self.velProxim[0]=filterdata(unfil[0],5)/1000
        self.velProxim[1]=filterdata(unfil[1],5)/1000
        return self.velProxim
    
    def accelProxim(self): #j1 is joint angle dataframe #n is number of datapoin
        BodyPart.velProxim(self)
        unfil=[np.zeros(self.sample),np.zeros(self.sample)]
        for i in range(1,self.sample-1):
            unfil[0][i] = (self.velProxim[0][i+1]-self.velProxim[0][i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
            unfil[1][i] = (self.velProxim[1][i+1]-self.velProxim[1][i-1])/(BodyPart.time[i+1]-BodyPart.time[i-1])
        self.accelProxim[0]=filterdata(unfil[0],5)
        self.accelProxim[1]=filterdata(unfil[1],5)
        return self.accelProxim
    
    def accelCoM(self):
        r=self.length-self.CoM
        BodyPart.accelProxim(self)
        BodyPart.alpha(self)
        unfil=[np.zeros(self.sample),np.zeros(self.sample)]
        for i in range(1,self.sample):
            unfil[0][i]=self.accelProxim[0][i]-(r*self.alpha[i]*math.sin(self.angle[i])+self.omega[i]*self.omega[i]*r*math.cos(self.angle[i]))
            unfil[1][i]=self.accelProxim[1][i]+r*self.alpha[i]*math.cos(self.angle[i])-self.omega[i]*self.omega[i]*r*math.sin(self.angle[i])
        self.accelCoM[0]=filterdata(unfil[0],5)
        self.accelCoM[1]=filterdata(unfil[1],5)
        return self.accelCoM    

    def Forces(self,R,M):
        BodyPart.accelCoM(self)
        for i in range(1,self.sample):
            self.ProximForce[0][i]=R[0][i]+self.mass*self.accelCoM[0][i];
            self.ProximForce[1][i]=R[1][i]+self.mass*self.accelCoM[1][i]+self.mass*9.81;
            self.ProximMoment[i]=M[i]+R[0][i]*self.length*math.sin(self.angle[i])-R[1][i]*self.length*math.cos(self.angle[i])-self.mass*9.81*self.CoM*math.cos(self.angle[i])+(self.Inertia+self.mass*self.CoM*self.CoM)*self.alpha[i];
        return  self.ProximForce,self.ProximMoment
    
def Power(omega1,omega2,M):
    Jw=omega2-omega1
    power=np.multiply(Jw,M)
    return power
    
def aCoM(m1,a1,m2,a2,m3,a3,m4,a4):
    acom=[np.zeros(106),np.zeros(106)]
    F=[np.zeros(106),np.zeros(106)]
    for i in range(0,2):
        for j in range(0,106):
            acom[i][j]=(m1*a1[i][j]+m2*a2[i][j]+m3*a3[i][j]+m4*a4[i][j])/(m1+m2+m3+m4)
            F[i][j]=(m1+m2+m3+m4)*(9.81+acom[i][j])
    return F

def filterdata(y,n):
    from scipy.signal import filtfilt
    b = [1.0 / n] * n
    a = 1
    yy = filtfilt(b,a,y)
    return yy;

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