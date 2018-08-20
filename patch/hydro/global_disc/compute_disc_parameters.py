import numpy as np
import scipy.integrate as integrate
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

def phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*2.0e33
    hb=p3*3.08e21
    Mds=p4*2.0e33
    Mdg=p5*2.0e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8

    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2+1e-10**2)*3.08e21
    
    phi_nfw=-4.0*np.pi*GG*rho_s*R_s**2*np.log(1+Rspher/R_s)/(Rspher/R_s)
    phi_bulge=-GG*Mb/np.sqrt(Rspher**2+hb**2)
    #phi_disc=-GG*(Mds+Mdg)/np.sqrt((rcirc+hr)**2+(np.abs(zcirc)+hz)**2)
    phi_disc=-GG*(Mds+Mdg)/np.sqrt(rcirc**2+(hr+np.sqrt(zcirc**2+hz**2))**2)
    phi_tot = phi_nfw + phi_bulge + phi_disc

    return phi_tot

def dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    dh = p6/1.0e4
    out = (phi(r+dh,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8))/(dh*3.08e21)

    return out
    
def dphi_dz(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hz=p7*3.08e21
    dh = p7/1e4
    out = (phi(r,z+dh,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8))/(dh*3.08e21)
    
    return out

def vcirc(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    dh = p7/1.0e4
    mp=1.66e-24
    kb=1.38062e-16

    rcirc = r*3.08e21
    
    dpdr=-2*kb*Tmu/mp/hr**2+dphi_dr(r,0.0,p0,p1,p2,p3,p4,p5,p6,p7,p8)
     
    out = (r*3.08e21*dpdr)**0.5

    return out
    
def kappa(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    dh = p6/1.0e4
    dpdr=dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    dpdr1=dphi_dr(r+dh,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    dpdr0=dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    term1 = 3*np.abs(dpdr)/(r*3.08e21)
    term2 = (np.abs(dpdr1)-np.abs(dpdr0))/(dh*3.08e21)
    out = (term1+term2)**0.5

    return out
    
def rho_scaling(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*2.0e33
    hb=p3*3.08e21
    Mds=p4*2.0e33
    Mdg=p5*2.0e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8
    mp=1.66e-24
    kb=1.38062e-16
    dh = p7/1.0e4
    
    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2)*3.08e21

    TT = Tmu
    
    out = np.exp(-(rcirc/hr)**2)*np.exp(-mp*(phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,0.0,p0,p1,p2,p3,p4,p5,p6,p7,p8))/kb/TT)
    
    return out

def Tmu_norm(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*2.0e33
    hb=p3*3.08e21
    Mds=p4*2.0e33
    Mdg=p5*2.0e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8
    mp=1.66e-24
    kb=1.38062e-16
    
    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2+1e-10**2)*3.08e21
    
    out = 0.99*mp/kb*hr**2/rcirc/2*dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
        
    return out

def total_integrand(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    out = 2*np.pi*rho_scaling(np.exp(r),np.exp(z),p0,p1,p2,p3,p4,p5,p6,p7,p8)*np.exp(2*r)*np.exp(z)
    return out

def vertical_integrand(z,r,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    out = rho_scaling(r,np.exp(z),p0,p1,p2,p3,p4,p5,p6,p7,p8)*np.exp(z)
    return out

#NFW parameters
h = 0.7 
H0 = h*100*1e5/3.08e24 #in cgs 
M200 = 1.5e11 # M200 in Msun 
R200 = (2*6.67e-8*2.0e33*M200/200/H0/H0)**(1.0/3.0)/3.08e21 # R200 in kpc
c200 = 8.03*(M200*h/1e12)**(-0.101) # Dutton and Maccio 2014 
R_s = R200/c200 # Scale radius in kpc
nfw_fac = np.log(1+c200)-c200/(1+c200)
rho_s = (M200*2.0e33)/(4.0*np.pi)/(R_s*3.08e21)**3/nfw_fac # Scale density in g/cm3
v200 = np.sqrt(6.67e-8*M200*2.0e33/R200/3.08e21)/1e5
T200 = 3.6e5*(v200/100)**2
print 'NFW Parameters: '
print 'H0 = ',h*100,' km/s/Mpc'
print 'M200 = ',M200,' Msun'
print 'R200 = ',R200,' kpc'
print 'c200 = ',c200
print 'T200 = ',T200
print 'R_s = ',R_s,' kpc'
print 'rho_s = ',rho_s,' g/cm3'
print ' '

#Bulge parameters
Mbulge = 1.0e-3  # Bulge mass in Msun
hbulge = 1.0e-6 # Bulge scale radius in kpc
print 'Bulge parameters: '
print 'Mbulge = ',Mbulge,' Msun'
print 'hbulge = ',hbulge,' kpc'
print ' '

#Disc parameters
Mdstar = 5.0e8 # Stellar disc mass in Msun
Mdgas = 1.0e9 # Stellar disc mass in Msun
hrdisc = 0.4 # Radial scale radius in kpc
hzdisc = 0.2 # Scale height in kpc
print 'Disc parameters: '
print 'Mdstar = ',Mdstar,' Msun'
print 'Mdgas = ',Mdgas,' Msun'
print 'hrdisc = ',hrdisc,' kpc'
print 'hzdisc = ',hzdisc,' kpc'
print ' '

#Compute gas disc normalization and central density
Tmu = 6.2e3
rmax=4.192*(3**0.33333)
Tmumax = Tmu_norm(rmax,0.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
params = (rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
normalization = integrate.nquad(total_integrand, [[np.log(1e-10*R200), np.log(2*R200)],[np.log(1.0e-10*R200), np.log(R200)]], args=params, opts = [{'limit' : 500}, {'limit' : 500}])
rho0gas = Mdgas*2.0e33/(2*normalization[0]*3.08e21**3)
print "rmax = ", rmax,' kpc'
print "rho0gas = ", rho0gas,' g/cm3'
print "nH0gas = ", 0.76*rho0gas/1.66e-24,' 1/cm3'
print "Tmax/mu = ",Tmumax," K"
print "T/mu = ",Tmu," K"
print ' '

#Define variables for plots
r = np.linspace(-3,np.log10(rmax),35)
z = np.linspace(-3,np.log10(rmax),35)
r = 10**r
z = 10**z
Sigmagas=np.zeros((len(r)))

#Compute Sigma_gas profile
i=0
for rr in r:
    params2 = (rr,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
    vals = integrate.quad(vertical_integrand,np.log(1.0e-12*R200), np.log(R200), args=params2, limit = 500)
    Sigmagas[i] = rho0gas*2*vals[0]*(1e3*3.08e18**3/2.0e33)+1e-20
    i=i+1

#Compute Qgas
dh = hzdisc/1e4
vc = vcirc(r,0.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)/1e5
c_s = (1.6666*1.38062e-16*Tmu/1.66e-24)**0.5
kap = kappa(r,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
GG = 6.67e-8
Qgas = kap*c_s/np.pi/GG/(Sigmagas*2.0e33/3.08e18**2)
#print Sigmagas, Qgas, kap

rhoa = rho_scaling(0.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)/rho_scaling(0.0,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
rhob = rho_scaling(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)/rho_scaling(1.0,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
rhoc = rho_scaling(2.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)/rho_scaling(2.0,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
rhor = rho_scaling(r,0.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)

drhodr = rho0gas*(rho_scaling(r+dh,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)-rho_scaling(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu))/(dh*3.08e21)
dPdr = 1.38e-16*drhodr*Tmu/1.66e-24
radial = -dPdr/rho0gas/rho_scaling(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu) -dphi_dr(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu) + (vc*1e5)**2/(r*3.08e21)
#print -dPdr/rho0gas/rho_scaling(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
#print -dphi_dr(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
#print (vc*1e5)**2/(r*3.08e21)
#print radial
#print -radial/dphi_dr(r,1.0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)

drhodz = rho0gas*(rho_scaling(1.0,z+dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)-rho_scaling(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu))/(dh*3.08e21)
dPdz = 1.38e-16*drhodz*Tmu/1.66e-24
vertical = -dPdz/rho0gas/rho_scaling(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu) -dphi_dz(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
#print -dPdz/rho0gas/rho_scaling(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
#print -dphi_dz(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
#print vertical
#print -vertical/dphi_dz(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)

plt.figure(0)
plt.axis([1e-3,1e1,1e-2, 1.0e4])
#plt.axis([1e-3,1e1,1e1, 1.0e10])
plt.xscale('log')
plt.yscale('log')
plt.plot(z,rhoa)
plt.plot(z,rhob)
plt.plot(z,rhoc)
plt.plot(r,rhor)
#plt.plot(r,Ta)
plt.plot(r,vc)
plt.plot(r,Sigmagas)
#plt.plot(r,Qgas)
#plt.plot(z,c/c*np.exp(-1))
plt.show()

