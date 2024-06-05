import matplotlib.pyplot as plt
import numpy
import random
from pandas import *
####DOubolel P recalcualte azeotrope

####Ternary plotting program
import matplotlib.pyplot as plt
import math
nipponnikal=math.sqrt(6)/2.
swin = math.sqrt(2)
def PointPlot(x1,x2,col):
    y = nipponnikal*(x1)
    x = (1-0.5*x1-x2 )*math.sqrt(2)
    plt.scatter(x,y,color = col)
    return 0
def LinePlot(x1,x2, lw = None, ck = None, ls = None):
    '''x1 and x2 are your points
        lw is your line width
        ck is your color
        ls is yoru linestyle
    '''
    if lw == None:
        lw = 1
    if ck == None:
        ck = 'b'
    if ls == None:
        ls = '-'
    x1 = ((1-0.5*x1[0]-x1[1] )*swin,nipponnikal*(x1[0]))
    x2 = ((1-0.5*x2[0]-x2[1] )*swin,nipponnikal*(x2[0]))
    plt.plot((x1[0],x2[0]),(x1[1],x2[1]), linewidth = lw, color= ck, linestyle = ls)
    return 0
def Scale():
    ###First part of scale constnat x2 lines
    x =[ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
        0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
        0.9 ,  0.95]
    for i in range(len(x)):
        LinePlot((x[i],1-x[i]),(0,1-x[i]), lw = 0.3, ls=':' )
        LinePlot((0,1-x[i]),(1-x[i],0),lw = 0.3, ls =':')
        LinePlot((x[i],1-x[i]),(x[i],0), lw=0.3, ls=':')
    return 0
def Show():
    width = 2
    plt.plot([0,math.sqrt(2)],[0,0],color = 'k',linewidth = width)##Bottom side
    plt.plot([0,math.sqrt(2)/2],[0,math.sqrt(6)/2.],color = 'k',linewidth = width)##left side
    plt.plot([math.sqrt(2),math.sqrt(2)/2],[0,math.sqrt(6)/2.],color = 'k', linewidth = width)
    ####Labeling
    plt.text(math.sqrt(2)/2-0.03,math.sqrt(6)/2.+0.03, '1', fontsize = 15)
    plt.text(-0.05,-0.03, '2', fontsize = 15)
    plt.text(math.sqrt(2)+0.02,-0.03, '3', fontsize = 15)
    plt.xticks([])
    plt.yticks([])
    plt.show()
    return 0




#########################
P= 101325
R=8.314


c12 = 5500
c23 = 5500
c13 = c31 = 4500


c12 = 3300
c23 = 2000
c13 = c31 = -2000

c12 = 5500
c23 = 5500
c13 = c31 = 4500




c12 = -4489.56
c23 = 3300.98
c13 = c31 = -4952.3



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
##Acetone-Choroform-Methanol System
v1 = 74.05
v2 = 80.67
v3 = 40.73


lam1211 = -349.3
lam2122 = -1586.4
lam1311 = -810.7
lam3133 = 2716.4
lam2322 = -1489.3
lam3233 = 7528.4


lam1211 = -349.3
lam2122 = -1586.4

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


    
def act1(x1,x2,x3,T):
    one = 1-numpy.log(x1+x2*A12(T)+x3*A13(T))
    two = (x1)/(x1+x2*A12(T)+x3*A13(T)) + (x2*A21(T))/(x1*A21(T)+x2+x3*A23(T)) + (x3*A31(T))/(x1*A31(T)+x2*A32(T)+x3)
    return numpy.exp(one-two)


def act2(x1,x2,x3,T):
    one = 1-numpy.log(x1*A21(T)+x2+x3*A23(T))
    two = (x1*A12(T))/(x1+x2*A12(T)+x3*A13(T))  + (x2)/(x1*A21(T)+x2+x3*A23(T))  +  (x3*A32(T))/(x1*A31(T)+x2*A32(T)+x3)
    return numpy.exp(one-two)


def act3(x1,x2,x3,T):
    one = 1-numpy.log(x1*A31(T)+x2*A32(T)+x3)
    two = (x1*A13(T))/(x1+x2*A12(T)+x3*A13(T))  + (x2*A23(T))/(x1*A21(T)+x2+x3*A23(T))  +  (x3)/(x1*A31(T)+x2*A32(T)+x3)
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

############
############
############
def dP1dT(T):
    return (P1(T+dT)-P1(T))/dT
def dP2dT(T):
    return (P2(T+dT)-P2(T))/dT
def dP3dT(T):
    return (P3(T+dT)-P3(T))/dT

#####################################
def lnact1(x1,x2,x3,T):
    return numpy.log(act1(x1,x2,x3,T))

def lnact2(x1,x2,x3,T):
    return numpy.log(act2(x1,x2,x3,T))

def lnact3(x1,x2,x3,T):
    return numpy.log(act3(x1,x2,x3,T))

def dlnact1dT(x1,x2,x3,T):
    return (lnact1(x1,x2,x3,T+dT)-lnact1(x1,x2,x3,T))/dT

def dlnact2dT(x1,x2,x3,T):
    return (lnact2(x1,x2,x3,T+dT)-lnact2(x1,x2,x3,T))/dT

def dlnact3dT(x1,x2,x3,T):
    return (lnact3(x1,x2,x3,T+dT)-lnact3(x1,x2,x3,T))/dT


#####################################

def Trial(x1,x2,oo):##maybe set Tguess to the average of the boiling points int he future
    x3 = oo#412 for the old exmpale
    Tguess = 56+273.15
    Told = 56+273.15
    #we are gonna try a newton raphson approach here instead of a try everything approach
    def G(T):
        return x1*act1(x1,x2,x3,T)*P1(T) + x2*act2(x1,x2,x3,T)*P2(T) + x3*act3(x1,x2,x3,T)*P3(T) - P
    def dGdT(T):
        return x1*act1(x1,x2,x3,T)*P1(T)*(dlnact1dT(x1,x2,x3,T)+(P1(T))**(-1)*dP1dT(T)) + x2*act2(x1,x2,x3,T)*P2(T)*(dlnact2dT(x1,x2,x3,T)+(P2(T))**(-1)*dP2dT(T)) + x3*act3(x1,x2,x3,T)*P3(T)*(dlnact3dT(x1,x2,x3,T)+(P3(T))**(-1)*dP3dT(T))
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
        print(x1,x2,x3)
        sqrt()
        return x1*(56+273.15)+x2*(31.2+273.15)+x3*(64.7+273.15)
    return Tguess


#####################################

####################################
def y1(x1,x2,x3):
        temp = Trial(x1,x2,x3)
        b = x1*act1(x1,x2,x3,temp)*P1(temp)/P
        if b>1.02 or b<0:
            print("y1 conc fail",x1,x2,x3,b)
            sqrt()
        return b
def y2(x1,x2,x3):
        temp = Trial(x1,x2,x3)
        b = x2*act2(x1,x2,x3,temp)*P2(temp)/P
        if b>1.02 or b<0:
            print("y2 conc fail",x1,x2,x3,b)
            sqrt()        
        return b
    
def y3(x1,x2,x3):########x3
        temp = Trial(x1,x2,x3)
        return x3*act3(x1,x2,x3,temp)*P3(temp)/P

#####################################
def dy1dx1(x1,x2,x3):
    return (y1(x1+dT,x2,x3-dT)-y1(x1,x2,x3))/dT
def dy1dx2(x1,x2,x3):
    return (y1(x1,x2+dT,x3-dT)-y1(x1,x2,x3))/dT
   
def dy2dx2(x1,x2,x3):
    return (y2(x1,x2+dT,x3-dT)-y2(x1,x2,x3))/dT

def dy2dx1(x1,x2,x3):
    return (y2(x1+dT,x2,x3-dT)-y2(x1,x2,x3))/dT

###################################################
def dis(PP1,PP2):
   # print(PP1)
  #  print(PP2)
    return numpy.sqrt((PP1[0]-PP2[0])**2+(PP1[1]-PP2[1])**2+(PP1[2]-PP2[2])**2)


def ZC(xstart):##Trying A(y-x)+f(c)
    dt = -0.1
    t=0
   # print(xstart)
    His = (xstart[0],xstart[1],xstart[2])
    while t>-200:
        ##Calculate A
        Ape = [[dy1dx1(xstart[0],xstart[1],xstart[2])-1,dy1dx2(xstart[0],xstart[1],xstart[2])],
               [dy2dx1(xstart[0],xstart[1],xstart[2]),dy2dx2(xstart[0],xstart[1],xstart[2])-1]]
        fc = [y1(xstart[0],xstart[1],xstart[2])-xstart[0] ,y2(xstart[0],xstart[1],xstart[2])-xstart[1]]


        dxdt = numpy.dot(Ape,fc)
        
        xstart = [dxdt[i]*dt+xstart[i] for i in range(len(dxdt))]

        xstart = [xstart[0],xstart[1],1-xstart[0]-xstart[1]]


        if numpy.isnan(xstart[0]):
            print("NAN Break")
            break
        if xstart[0] > 1 or xstart[0] < 0 or xstart[1] > 1 or xstart[1] < 0 or xstart[2]>1 or xstart[2]<0:
            print("conc fails")
            break
        if xstart[0]+xstart[1]+xstart[2]>1.05:
            print("sumfails")
            break
        if dis(His,xstart)<1e-5:
            print("converge")
            break
       # LinePlot((xstart[0],xstart[1]),(His[0],His[1]),ck=c)
        His = (xstart[0],xstart[1],xstart[2])
        t=dt+t
       # print (dTdx1(xstart[0],xstart[1])*(y1(xstart[0],xstart[1])-xstart[0]), dTdx2(xstart[0],xstart[1])*(y2(xstart[0],xstart[1])-xstart[1]))
  #  PointPlot(xstart[0],xstart[1],"b")
    return (xstart[0],xstart[1],xstart[2])



Aze = (0.34005193207098233, 0.22167654315899626, 1-0.22167654315899626-0.34005193207098233)

def stabFinder():
    x = numpy.linspace(0,1,100)
    y = numpy.linspace(0,1,100)


    for i in x:
        for j in y:
            x3 = 1-i-j
            if x3<0:
                continue
           # PointPlot(i,j,'g')
            if dis(ZC((i,j,1-i-j)),Aze)<1e-1:
                PointPlot(i,j,"y")

    

    Scale()

    return 0

stabFinder()




PointPlot(Aze[0],Aze[1],"r")

plt.xlim(-0.1,1.45)
plt.ylim(-0.1,1.45)

plt.title("Acetone-Chloroform-Methanol Ternary System, 1atm, Area of Stability")

Scale()
Show()




def azefinder():
    x = numpy.linspace(0,1,100)
    y = numpy.linspace(0,1,100)

    for i in x:
        for j in y:
            x3 = 1-i-j
            if x3<0:
                continue
            y11 = y1(i,j,x3)
            y22 = y2(i,j,x3)
            y33=1-y11-y22
            if dis((i,j,x3),(y11,y22,y33))<1e-3:
                PointPlot(i,j,"g")

    PointPlot(0.381,0.241,"r")

    Scale()
    Show()

    return 0



