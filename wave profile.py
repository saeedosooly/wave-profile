import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
###
##############################
######## Linear wave theory
##################################
t = 8
x = np.linspace(0,10,101)
h = 0.60
H = 0.08
print(t)
a = H/2 #m
K = (2*math.pi)/6.95
y = a*(np.cos(K*x - (2.09*t)))

################################################################3
##second - order Stokes wave
################################################
y2_1 = a*(np.cos(K*x - (2.09*t)))
y2_2_1 = (H*H*K/16)*((np.cosh(K*h))/(math.pow(np.sinh(K*h),3)))
y2_2_2 = ((2 + np.cosh(2*K*h))) * (np.cos(2*(K*x - (2.09*t))))

y2_2 = y2_2_1 * y2_2_2
y2 = y2_1 + y2_2

#######################################
##Third-order Stokes wave 
##################

y3_1= ((1 -((K*a)**2)/16))*np.cos(K*x - (2.09*t))
y3_2 = 0.5*K*a*np.cos(2*(K*x - (2.09*t)))
y3_3 = ((3/8)*((K*a)**2))*np.cos(3*(K*x - (2.09*t)))
y3 = a*(y3_1 +y3_2 +y3_3)

######################################################
######### cnoidal wave#######Korteweg–de Vries equation method
##################################
m = 1 - (math.pow(10,-0.5)) #
K__m = special.ellipk(m) #Complete elliptic integral of the first kind
E__m = special.ellipe(m) #Complete elliptic integral of the second kind
c = (math.sqrt(9.81*h))*(1 + ((H/(m*h))*(1 - (m/2) - 1.5*(E__m/K__m))))
landa = h*K__m*math.sqrt((16/3)*(m*h)/H)
delta = landa/(2*K__m)
eta2 = (H/m)*(1 - m - (E__m/K__m))
#print(y,y2,y3)
j = special.ellipj((x - (c*t))/delta,m) #Jacobian elliptic functions
cn = j[1]
eta = (cn*cn*H) + eta2   
####################
####### Graph #######
#u = [0 for x in range(1000)]
plt.plot(x,y+.6,x,y2+.6,x,y3+.6,x,eta+.6)
plt.legend(["Linear wave theory","second-order Stokes wave", "Third-order Stokes wave","cnoidal wave"],loc='best')
plt.xlabel(' x(m)')
plt.ylabel('\u03B7 (m)')
plt.grid()
plt.ylim(0,.7)
plt.title('Linear,second-order,Third-order Stokes,and cnoidal wave profile in x = 0-10 m & t= 5-8 s',loc = 'left')

#pylab.get_current_fig_manager().window.showMaximized()
#figManager = plt.get_current_fig_manager()
#figManager.resize(*figManager.window.maxsize())
#plt.savefig('x=0-30 m,t=5-8 s.png')
plt.show()
