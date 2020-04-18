# generate random integer values
from random import seed, random
from random import randint
import cv2
import numpy as np
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets,QtCore,QtGui,uic 
import pylab 
import time
import os
# print(os.path.isfile('c_1.png'))
img = cv2.imread('c_1.png',0)
xdata = []
ydata = []  
zdata =[]
tdata=[]
dT = 1	
T = 500
N = int(np.ceil(T/dT)+1)
df = 10	
T1 = 600	
T2 = 100


def fourier():
    f = np.fft.fft2(img)
    fshift = np.fft.fftshift(f)
    magnitude_spectrum =20*np.log(np.abs(fshift))
    plt.subplot(121)
    plt.imshow(img, cmap='gray')
    plt.title('Input Image'), plt.xticks([]), plt.yticks([])
    plt.subplot(122)
    plt.imshow(magnitude_spectrum, cmap='gray')
    plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
    plt.show()
    return


def randomNums():
    seed(1)
    z = np.arange(0, 2000, 1) 
    B=[]
    for _ in range(0,2000):

        B.append(random()*0.5+3)
    axes = pylab.gca()
    axes.set_xlim(0,2000)
    axes.set_ylim(-3,5)
    plt.plot(z, B,label='B')
    plt.xlabel('Z in mm')
    plt.ylabel('B in Tesla')
    plt.legend()
    plt.show()
    return


def zrot(phi):
   
    Rz = [[np.cos(phi) ,-np.sin(phi), 0],[np.sin(phi) ,np.cos(phi), 0],[ 0, 0 ,1]]
    return(Rz)

def freeprecess(dT ,T1 ,T2 , df):
    phi = 2*np.pi*df*dT/1000
    E1 = np.exp(-dT/T1)	
    E2 = np.exp(-dT/T2)
 
    Afp = np.dot([[E2,0,0],[0,E2,0],[0,0,E1]],zrot(phi))
      
    Bfp = [0 ,0 ,1-E1]
   
    return(Afp,Bfp)
    
def plot(axes,DataX,DataY,DataZ,timeData):
    ydata=[]
    xdata=[]
    if (DataZ == []):
        line, = axes.plot(xdata, ydata, 'r-')

       
        axes.set_xlim(min(DataX)+0.5*min(DataX),max(DataX)+0.1*max(DataX))
        axes.set_ylim(min(DataY)+0.5*min(DataY),max(DataY)+0.1*max(DataY))
        pylab.xlabel('MX')
        pylab.ylabel('MY')

        for i in range(N):
            xdata.append(DataX[i])
            ydata.append(DataY[i])

            line.set_xdata(xdata)

            line.set_ydata(ydata)
            plt.draw()
            plt.pause(1e-17)
            time.sleep(0.001)
        plt.show()
    else:
        line, = axes.plot(xdata, ydata, 'r-')
        
        axes.set_xlim(0,1000)
        axes.set_ylim(-0.1,1)  
     
        pylab.xlabel('Time')
        pylab.ylabel('Magnetization')

        for i in range(1,N):
            xdata.append(DataX[i-1])
            ydata.append(DataY[i-1])  
            zdata.append(DataZ[i-1])
            pylab.plot(timeData[:i],xdata,color = [0.9290, 0.6940, 0.1250],linewidth=1.5)
            pylab.plot(timeData[:i],ydata,color = [0.8500, 0.3250, 0.0980],linewidth=1.5)
            pylab.plot(timeData[:i],zdata,color = 'purple',linewidth=1.5)
            pylab.plot(timeData[:i],np.exp(-timeData[:i]/T2),color =[0.25, 0.25, 0.25],linewidth=1.5) 
            pylab.draw()
            pylab.pause(1e-117)
        
        plt.show()
    return



def blochEquation():

    
    A,B = freeprecess(dT,T1,T2,df)
 
    M = np.empty((N,3))
    
    M[0,:] =np.array([1,0,0])
    

    for k in range (N-1):

        M[k+1,:] = np.dot(A,M[k,:])+B
  

    pylab.subplot(111)

    DataX= M[:,0]

    DataY=M[:,1]
    DataZ = M[:,2]
    timeData = dT * np.arange(N)
    axes = pylab.gca()
   
    plot(axes,DataX,DataY,[],[])
    plot(axes,DataX,DataY,DataZ,timeData)










fourier()
randomNums()
blochEquation()
