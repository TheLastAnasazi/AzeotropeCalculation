import matplotlib.pyplot as plt
import numpy
import random
from pandas import *
####DOubolel P recalcualte azeotrope


import math
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
 
fig = plt.figure()
ax = plt.axes(projection ='3d')
 



#########################
P= 101325
R=8.314





c12 = 3300
c23 = 2000
c13 = c31 = -2000

c12 = 5500
c23 = 5500
c13 = c31 = 4500




c12 = -4489.56
c23 = 3300.98
c13 = c31 = -4952.3

c12 = 5500
c23 = 5500
c13 = c31 = 4500



dT = 0.000001

##############################


    ####1 con Margules
def act1(x1,x2,x3,T):
        one = c12*(x2**2+x2*x3)
        two = -1*c23*x2*x3
        three = c13*(x3**2+x2*x3)
        return numpy.exp((R*T)**(-1)*(one + two+three))

def act2(x1,x2,x3,T):
        one = c23*(x3**2+x3*x1)
        two = c12*(x1*x3+x1**2)
        three = -1*c13*(x3*x1)
        return numpy.exp((R*T)**(-1)*(one + two+three))

def act3(x1,x2,x3,T):
        one = -1*c12*x1*x2
        two = c23*(x1*x2+x2**2)
        three = c13*(x2*x1+x1**2)
        return numpy.exp((R*T)**(-1)*(one + two+three))



    

#### Wilson Model
##Acetone-Choroform-Methanol-cyclohexane System
v1 = 74.05##acetone
v2 = 80.67##chlrofomr
v3 = 40.73##methanol
v4 = 108.7#cyclehexane


lam1211 = -349.3##acetone-chloro
lam2122 = -1586.4
lam1311 = -810.7##Acetone methanol
lam3133 = 2716.4
lam2322 = -1489.3###Chloro methanol
lam3233 = 7528.4
lam1411 = 2795.99556532####Acetone cyclhexane
lam4144 = 2691.855771214645


lam2422 = 1##chloro cyclohexane 2-4 willbe treated as ideal.Need to figure out what that means still lol
lam4244 = 1##ideal case is they are both 1 which is interestint its not0

lam3433 = 4895.391728008291###methanol cyclehexane
lam4344 = 8605.952012184827

#lam1311 = 2795.99556532
#lam3133 = 2691.855771214645
#lam2322 = 1
#lam3233 = 1


#lam1211 = 2795.99556532
#lam2122 = 2691.855771214645

def A12(T):
    one = v2/v1
    return one*numpy.exp(-1*(lam1211)/(R*T))

def A21(T):
    one = v1/v2
    return one*numpy.exp(-1*(lam2122)/(R*T))

def A13(T):
    one = v3/v1
    return one*numpy.exp(-1*(lam1311)/(R*T))

def A31(T):
    one = v1/v3
    return one*numpy.exp(-1*(lam3133)/(R*T))

def A23(T):
    one = v3/v2
    return one*numpy.exp(-1*(lam2322)/(R*T))

def A32(T):
    one = v2/v3
    return one*numpy.exp(-1*(lam3233)/(R*T))
#####adding comp 4
def A14(T):
    one = v4/v1
    return one*numpy.exp(-1*(lam1411)/(R*T))

def A41(T):
    one = v1/v4
    return one*numpy.exp(-1*(lam4144)/(R*T))

def A24(T):
    one = v4/v2
    return one*numpy.exp(-1*(lam2422)/(R*T))

def A42(T):
    one = v2/v4
    return one*numpy.exp(-1*(lam4244)/(R*T))

def A34(T):
    one = v4/v3
    return one*numpy.exp(-1*(lam3433)/(R*T))

def A43(T):
    one = v3/v4
    return one*numpy.exp(-1*(lam4344)/(R*T))



    
def act1(x1,x2,x3,x4,T):
    D1 = x1 +        x2*A12(T) + x3*A13(T) + x4*A14(T)
    D2 = x1*A21(T) + x2        + x3*A23(T) + x4*A24(T)
    D3 = x1*A31(T) + x2*A32(T) + x3        + x4*A34(T)
    D4 = x1*A41(T) + x2*A42(T) + x3*A43(T) + x4
    one = 1-numpy.log(x1+x2*A12(T)+x3*A13(T)+x4*A14(T))
    two = (x1)/(D1) + (x2*A21(T))/(D2) + (x3*A31(T))/(D3)+(x4*A41(T))/D4
    return numpy.exp(one-two)


def act2(x1,x2,x3,x4,T):
    D1 = x1 +        x2*A12(T) + x3*A13(T) + x4*A14(T)
    D2 = x1*A21(T) + x2        + x3*A23(T) + x4*A24(T)
    D3 = x1*A31(T) + x2*A32(T) + x3        + x4*A34(T)
    D4 = x1*A41(T) + x2*A42(T) + x3*A43(T) + x4
    one = 1-numpy.log(x1*A21(T)+x2+x3*A23(T)+x4*A24(T))
    two = (x1*A12(T))/(D1)  + (x2)/(D2)  +  (x3*A32(T))/(D3)+ (x4*A42(T))/D4
    return numpy.exp(one-two)


def act3(x1,x2,x3,x4,T):
    D1 = x1 +        x2*A12(T) + x3*A13(T) + x4*A14(T)
    D2 = x1*A21(T) + x2        + x3*A23(T) + x4*A24(T)
    D3 = x1*A31(T) + x2*A32(T) + x3        + x4*A34(T)
    D4 = x1*A41(T) + x2*A42(T) + x3*A43(T) + x4
    one = 1-numpy.log(x1*A31(T)+x2*A32(T)+x3+x4*A34(T))
    two = (x1*A13(T))/(D1)  + (x2*A23(T))/(D2)  +  (x3)/(D3)+ (x4*A43(T))/D4
    return numpy.exp(one-two)

def act4(x1,x2,x3,x4,T):
    D1 = x1 +        x2*A12(T) + x3*A13(T) + x4*A14(T)
    D2 = x1*A21(T) + x2        + x3*A23(T) + x4*A24(T)
    D3 = x1*A31(T) + x2*A32(T) + x3        + x4*A34(T)
    D4 = x1*A41(T) + x2*A42(T) + x3*A43(T) + x4
    one = 1-numpy.log(x1*A41(T)+x2*A42(T)+x3*A43(T)+x4)
    two = (x1*A14(T))/(D1)  + (x2*A24(T))/(D2)  +  (x3*A34(T))/(D3)+ (x4)/D4
    return numpy.exp(one-two)




########################################

########################################

def gen(C,T):
        return numpy.exp(C[0]+C[1]/T + C[2]*numpy.log(T) + C[3]*T**C[4])

def P1(T):#IPA
        return gen([76.964,-7623.8,-7.4924,5.9436*10**(-18),6],T)

def P2(T):#Water
        return gen([73.649,-7258.2,-7.3037,4.1653*10**(-6),2],T)

def P3(T):#2butanol
        return gen([152.54,-11111,-19.025,1.0426*10**(-5),2],T)






###################

def gen(C,T):## log(Psat) = A-B/(T+C)###T in Celcius##Psat in torr
        return 133.322*10**(C[0]-C[1]/(T+C[2]-273.15))
    
def P1(T):#Acetone
        return gen([7.11714,1210.595,229.664],T)

def P2(T):#Chloroform
        return gen([6.95465,1170.966,226.232],T)

def P3(T):#Methanol
        return gen([8.08097,1582.271,239.726],T)

def P4(T):#Cyclohexane##https://webbook.nist.gov/cgi/inchi?ID=C110827&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
    return 100000/133.322*gen([4.13983,1316.554,-35.581+273.15],T)


############
############
############
def dP1dT(T):
    return (P1(T+dT)-P1(T))/dT
def dP2dT(T):
    return (P2(T+dT)-P2(T))/dT
def dP3dT(T):
    return (P3(T+dT)-P3(T))/dT
def dP4dT(T):
    return (P3(T+dT)-P3(T))/dT

#####################################
def lnact1(x1,x2,x3,x4,T):
    return numpy.log(act1(x1,x2,x3,x4,T))

def lnact2(x1,x2,x3,x4,T):
    return numpy.log(act2(x1,x2,x3,x4,T))

def lnact3(x1,x2,x3,x4,T):
    return numpy.log(act3(x1,x2,x3,x4,T))

def lnact4(x1,x2,x3,x4,T):
    return numpy.log(act4(x1,x2,x3,x4,T))

def dlnact1dT(x1,x2,x3,x4,T):
    return (lnact1(x1,x2,x3,x4,T+dT)-lnact1(x1,x2,x3,x4,T))/dT

def dlnact2dT(x1,x2,x3,x4,T):
    return (lnact2(x1,x2,x3,x4,T+dT)-lnact2(x1,x2,x3,x4,T))/dT

def dlnact3dT(x1,x2,x3,x4,T):
    return (lnact3(x1,x2,x3,x4,T+dT)-lnact3(x1,x2,x3,x4,T))/dT

def dlnact4dT(x1,x2,x3,x4,T):
    return (lnact4(x1,x2,x3,x4,T+dT)-lnact4(x1,x2,x3,x4,T))/dT


#####################################

def Trial(x1,x2,x3,x4):##maybe set Tguess to the average of the boiling points int he future
    #x3 = oo#412 for the old exmpale
    Tguess = 340#56+273.15
    Told = 402
    #we are gonna try a newton raphson approach here instead of a try everything approach
    def G(T):
        return x1*act1(x1,x2,x3,x4,T)*P1(T) + x2*act2(x1,x2,x3,x4,T)*P2(T) + x3*act3(x1,x2,x3,x4,T)*P3(T) + x4*act4(x1,x2,x3,x4,T)*P4(T) - P
    def dGdT(T):
        return x1*act1(x1,x2,x3,x4,T)*P1(T)*(dlnact1dT(x1,x2,x3,x4,T)+(P1(T))**(-1)*dP1dT(T))   +    x2*act2(x1,x2,x3,x4,T)*P2(T)*(dlnact2dT(x1,x2,x3,x4,T)+(P2(T))**(-1)*dP2dT(T))     +    x3*act3(x1,x2,x3,x4,T)*P3(T)*(dlnact3dT(x1,x2,x3,x4,T)+(P3(T))**(-1)*dP3dT(T))   +    x4*act4(x1,x2,x3,x4,T)*P4(T)*(dlnact4dT(x1,x2,x3,x4,T)+(P4(T))**(-1)*dP4dT(T))
    for i in range(20):
        #print Tguess,i
        Tguess -= G(Tguess)/dGdT(Tguess)
        if Told-Tguess<1e-4:
            break
        else:
            Told = Tguess
        #print(Tguess)
      #  print(x1,x2,x3,Tguess)
    if numpy.isnan(Tguess):
        print(x1,x2,x3,x4)
        return x1*(56+273.15)+x2*(31.2+273.15)+x3*(64.7+273.15)
    return Tguess


#####################################

####################################
def y1(x1,x2,x3,x4):
        temp = Trial(x1,x2,x3,x4)
        b = x1*act1(x1,x2,x3,x4,temp)*P1(temp)/P
        if b>1.001 or b<0:
            print("y1 conc fail",x1,x2,x3,x4,temp)
            sqrt()
        return b
def y2(x1,x2,x3,x4):
        temp = Trial(x1,x2,x3,x4)
        b = x2*act2(x1,x2,x3,x4,temp)*P2(temp)/P
        if b>1.001 or b<0:
            print("y2 conc fail",x1,x2,x3,x4,temp)
            sqrt()        
        return b
    
def y3(x1,x2,x3,x4):########x3
        temp = Trial(x1,x2,x3,x4)
        return x3*act3(x1,x2,x3,x4,temp)*P3(temp)/P

def y4(x1,x2,x3,x4):########x3
        temp = Trial(x1,x2,x3,x4)
        return x4*act4(x1,x2,x3,x4,temp)*P4(temp)/P



    

#####################################
def dy1dx1(x1,x2,x3,x4):
    return (y1(x1+dT,x2,x3,x4-dT)-y1(x1,x2,x3,x4))/dT
def dy1dx2(x1,x2,x3,x4):
    return (y1(x1,x2+dT,x3,x4-dT)-y1(x1,x2,x3,x4))/dT
def dy1dx3(x1,x2,x3,x4):
    return (y1(x1,x2,x3+dT,x4-dT)-y1(x1,x2,x3,x4))/dT


  
def dy2dx2(x1,x2,x3,x4):
    return (y2(x1,x2+dT,x3,x4-dT)-y2(x1,x2,x3,x4))/dT
def dy2dx1(x1,x2,x3,x4):
    return (y2(x1+dT,x2,x3,x4-dT)-y2(x1,x2,x3,x4))/dT
def dy2dx3(x1,x2,x3,x4):
    return (y2(x1,x2,x3+dT,x4-dT)-y2(x1,x2,x3,x4))/dT

def dy3dx1(x1,x2,x3,x4):
    return (y3(x1+dT,x2,x3,x4-dT)-y3(x1,x2,x3,x4))/dT
def dy3dx2(x1,x2,x3,x4):
    return (y3(x1,x2+dT,x3,x4-dT)-y3(x1,x2,x3,x4))/dT
def dy3dx3(x1,x2,x3,x4):
    return (y3(x1,x2,x3+dT,x4-dT)-y3(x1,x2,x3,x4))/dT




###################################################
def dis(PP1,PP2):
   # print(PP1)
  #  print(PP2)
    return numpy.sqrt((PP1[0]-PP2[0])**2+(PP1[1]-PP2[1])**2+(PP1[2]-PP2[2])**2)
###################################################


def fc(xstart,c):
    #if 1==1:
   #     return 0
    origin  =xstart
    dt = -0.1
    t=0
    xstart = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]
    His = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]
    while t<40:
        wi = y1(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        wo = y2(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        we =  y3(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        dy1 = (wi-xstart[0])*dt
        dy2 = (wo-xstart[1])*dt
        dy3 = (we-xstart[2])*dt
      #  PointPlot(xstart[0],xstart[1],c)
      #  print xstart, wi, y1lin(xstart[0],xstart[1])
        xstart[0] = xstart[0] + dy1
        xstart[1] = xstart[1] + dy2
        xstart[2] = xstart[2] + dy3
        xstart[3] = 1-xstart[0] - xstart[1] - xstart[2]
        
   #     if dis(His,(xstart[0], xstart[1], 1-xstart[0]-xstart[1]))<1e-4:
    #        break
        ax.plot3D([xstart[0],His[0]],[xstart[1],His[1]],[xstart[2],His[2]],c=c)
      #  LinePlot((xstart[0],xstart[1]),(His[0],His[1]),ck=c)
        His = (xstart[0],xstart[1],xstart[2],xstart[3])
      #  print xstart, His #, dy1
        t-=dt
    t=0
    dt = 0.1
    xstart= origin
    xstart = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]
    His = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]
    while t<40:
        wi = y1(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        wo = y2(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        we =  y3(xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2])
        dy1 = (wi-xstart[0])*dt
        dy2 = (wo-xstart[1])*dt
        dy3 = (we-xstart[2])*dt
      #  PointPlot(xstart[0],xstart[1],c)
      #  print xstart, wi, y1lin(xstart[0],xstart[1])
        xstart[0] = xstart[0] + dy1
        xstart[1] = xstart[1] + dy2
        xstart[2] = xstart[2] + dy3
        xstart[3] = 1-xstart[0] - xstart[1] - xstart[2]
        
   #     if dis(His,(xstart[0], xstart[1], 1-xstart[0]-xstart[1]))<1e-4:
    #        break
        ax.plot3D([xstart[0],His[0]],[xstart[1],His[1]],[xstart[2],His[2]],c=c)
      #  LinePlot((xstart[0],xstart[1]),(His[0],His[1]),ck=c)
        His = (xstart[0],xstart[1],xstart[2],xstart[3])
      #  print xstart, His #, dy1
        t=t+dt
    #ax.scatter(His[0],His[1],His[2],c="b")
    return 0#xstart[0], xstart[1]


###########################################################



def ZC(xstart,c):##Trying A(y-x)+f(c)
    dt = -0.1
    t=0
    print(xstart)
    xstart = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]
    His = (xstart[0],xstart[1],xstart[2],xstart[3])
    while t>-200:
        ##Calculate A
        Ape = [[dy1dx1(xstart[0],xstart[1],xstart[2],xstart[3])-1,dy1dx2(xstart[0],xstart[1],xstart[2],xstart[3]),dy1dx3(xstart[0],xstart[1],xstart[2],xstart[3])],
               
               [dy2dx1(xstart[0],xstart[1],xstart[2],xstart[3]),dy2dx2(xstart[0],xstart[1],xstart[2],xstart[3])-1,dy2dx3(xstart[0],xstart[1],xstart[2],xstart[3])],

               [dy3dx1(xstart[0],xstart[1],xstart[2],xstart[3]),dy3dx2(xstart[0],xstart[1],xstart[2],xstart[3]),dy3dx3(xstart[0],xstart[1],xstart[2],xstart[3])-1]
               ]
        fc = [y1(xstart[0],xstart[1],xstart[2],xstart[3])-xstart[0] ,y2(xstart[0],xstart[1],xstart[2],xstart[3])-xstart[1], y3(xstart[0],xstart[1],xstart[2],xstart[3])-xstart[2]]


        dxdt = numpy.dot(Ape,fc)
        
        xstart = [dxdt[i]*dt+xstart[i] for i in range(len(dxdt))]

        xstart = [xstart[0],xstart[1],xstart[2],1-xstart[0]-xstart[1]-xstart[2]]


        if numpy.isnan(xstart[0]):
            print("NAN Break")
            break
        if xstart[0] > 1 or xstart[0] < 0 or xstart[1] > 1 or xstart[1] < 0 or xstart[2]>1 or xstart[2]<0:
            print("conc fails")
            break
        if xstart[0]+xstart[1]+xstart[2]>1.05:
            print("sumfails")
            break
        if dis(His,xstart)<1e-7:
            #PointPlot(xstart[0],xstart[1],"b")
            ax.scatter(His[0],His[1],His[2],c="b")
            print("converge")
            break
        ax.plot3D([xstart[0],His[0]],[xstart[1],His[1]],[xstart[2],His[2]],c="k")
       # LinePlot((xstart[0],xstart[1]),(His[0],His[1]),ck=c)
        His = (xstart[0],xstart[1],xstart[2],xstart[3])
        t=dt+t
       # print (dTdx1(xstart[0],xstart[1])*(y1(xstart[0],xstart[1])-xstart[0]), dTdx2(xstart[0],xstart[1])*(y2(xstart[0],xstart[1])-xstart[1]))
    
    #print("Done")
    
    return xstart[0],xstart[1],xstart[2],xstart[3]

###################











def ksh(y,cxc):
  #  ax.scatter(y[0],y[1],y[2],c='g')
    aa = fc([y[0],y[1],y[2]],cxc)
    bb = ZC([y[0],y[1],y[2]],cxc)
    print ('thre3e',aa,bb)
    return 0



x1aze = 0.39
x2aze = 0.32

colist = ["m","b",'r','g','y','c']
import random
for i in range(15):
    a = random.random()
    b = random.random()
    c = random.random()
    d = random.random()
    norm = a+b+c+d
    a = a/norm
    b = b/norm
    c = c/norm
    d = d/norm
    ksh((a,b,c),colist[random.randint(0,len(colist)-1)])



ax.plot3D([0,0,1,0],[0,1,0,0],[1,0,0,1],"k")
ax.scatter(0,0,0,c="r")
ax.scatter(1,0,0,c="r")
ax.scatter(0,1,0,c="r")
ax.scatter(0,0,1,c="r")
ax.set_xlabel('Acetone', fontsize=12)
ax.set_ylabel('Chloroform', fontsize=12)
ax.set_zlabel('Methanol', fontsize=12)

ax.set_title("Acetone Chloroform Methanol Cyclohexane System, 1atm")

plt.show()


ksh((0.39,0.32,0.29),"m")

#sqrt()
ksh((0.1,0.8,0.1),'b')
ksh((0.5,0.4,0.1),"r")
ksh((0.05,0.15,0.8),'r')
ksh((0.05,0.9,0.05),"y")
ksh((0.9,0.05,0.05),"r")
ksh((0.20,0.50,0.3),"g")

ksh((0.2,0.2,0.6),"r")
ksh((0.3,0.3,0.4),"b")
ksh((0.5,0.1,0.4),"y")
ksh((0.4,0.1,0.5),"m")
ksh((0.381,0.241,1-0.381-0.241),"g")
ksh((0.331,0.201,1-0.331,0.201),"b")
ksh((0.3,0.25,0.45),"y")
ksh((0.75,0.1,0.15),"r")
ksh((0.65,0.15,0.2),"b")
ksh((0.55,0.25,0.2),"b")



Scale()
Show()









##########################################################

def binTxy():
    x = numpy.linspace(0,1,1000)
    yy = [i*act1(i,1-i,0,Trial(i,1-i,0))*P1(Trial(i,1-i,0))/P for i in x]
    TT = [Trial(i,1-i,0) for i in x]

    plt.plot([0,1],[0,1])
    plt.scatter(x,yy,s=1,color='r')
    plt.show()

    plt.scatter(x,TT,s=1)
    plt.scatter(yy,TT,s=1)
    plt.show()
    return 0



x1aze = 0.605
x2aze = 1-0.605
TemAze = 327.74

Alpha1 = numpy.log(P/P3(TemAze))
Alpha2 = numpy.log(P/P4(TemAze))


#from scipy.optimize import fsolve
#import math

#def equations(p):
#    A12, A21 = p
#3    Volk = (A12)/(x1aze+A12*x2aze)-(A21)/(x2aze+A21*x1aze)
#3    return (-numpy.log(x1aze+A12*x2aze)+x2aze*(Volk)-Alpha1, -numpy.log(x2aze+A21*x1aze)-x1aze*Volk - Alpha2)

#x, y =  fsolve(equations, (0.5, 0.5))

#print(equations((x, y)))

#A14=0.24395213811337438
#A41 = 0.5461664149360087

#xx = numpy.linspace(0,1,100)
#yy = numpy.linspace(0,1,100)
#Piple = [[],[]]
#for i in xx:
#    for j in yy:
#        ff = equations([i,j])
#        Piple[0].append(ff[0])
#        Piple[1].append(ff[1])
#    plt.scatter(Piple[0],Piple[1])
#    Piple = [[],[]]
#plt.scatter(Piple[0],Piple[1])
#plt.scatter(0,0,color="r")
#plt.show()







