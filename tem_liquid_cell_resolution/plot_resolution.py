import numpy as np
import matplotlib.pyplot as plt

c = 3*10**8 # speed of light (m/s)
h = 6.626*10**(-34) # planck's constant (kg m^2 s^-1)
E_0 = 511 # electron rest energy (keV)
N_A = 6.022*10**23 # Avagadro's number
e = 1.602*10**(-19) # charge of electron (C)
a_h = 5.29*(10**(-11)) # bohr radius (m)
epsilon_0 = 8.854*10**(-15) # permittivity of free space (C^2 g^-1 m^-3 s^2)

E = 300 # incident electron energy (keV)
Lambda = 1.9687*10**(-12) # relativistic wavelength of electron (m) https://www.jeol.com/words/emterms/20121023.071258.php#gsc.tab=0
C_c = 2.4 # chromatic aberration coefficient (mm)
C_s = 1.6 # spherical aberration coefficient (mm)
d_objective_aperture = 60*10**(-6) # diameter of objective aperture (m)
alpha_0 = np.arcsin(((Lambda)/(2*d_objective_aperture))) # half-angle of objective aperture

v = c*(1-(1/((1+(E/E_0))**2)))**(1/2) # relativistic velocity (m/s)

H = np.arange(0*10**(-7), 2000*10**(-7), 50*10**(-7)) #liquid thickness (cm) #deltaE should have H in cm, but why?

t_SiN = 20*10**(-9) # membrane thickness (m)

R_O = 4.8*(10**(-11)) # radius of oxygen atom (m)
R_H = 5.3*(10**(-11)) # radius of hydrogen atom (m)
R_Si = 11.1*(10**(-11)) # radius of silicon atom (m)
R_N = 5.6*(10**(-11)) # radius of nitrogen atom (m)
#R_C = 6.7*(10**(-11)) # radius of carbon atom (m)
#R_Zr = 20.6*(10**(-11)) # radius of zirconium atom (m)

flux = 0.1 # dose rate (electrons/Angstrom^2/s)
exposure_time = 1 # exposure time (s)
probe_diameter = 12.59*10**4 # beam/probe diameter (Angstrom)
probe_radius = probe_diameter/2 # beam/probe radius (Angstrom)
probe_area = np.pi*((probe_radius)**2) # beam/probe area (Angstrom^2)
N_0 = flux*exposure_time*probe_area # number of incident electrons 


def sigma(x, y, R_0, R_1):
    sigma = (x/(x+y))*(np.pi*(2*R_0)**2) + (y/(y+x))*(np.pi*(2*R_1)**2)
    return sigma

def l(A, sigma, rho):
    l = A/(sigma*rho*N_A)
    return l


A_water = ((2*1.008)+(1*15.999))*1.674*10**(-24)# atomic weight of water (g)
A_SiN = ((3*28.084)+(4*14.0064))*1.674*10**(-24) # atomic weight of SiN (g)
rho_object = 1.237*10**6 # specimen density of UiO-66 (g/m^3)
rho_water = 1.0*10**6 # specimen density of water (g/m^3)
rho_SiN = 3.17*10**6 # specimen density of SiN (g/m^3)

Z_gold = 79 # atomic number of Gold
A_gold = 196.967 # atomic weight of Gold (g)
R_gold = 17.4*10**(-11) # radius of gold atom (m)


l_o = l(A_gold, sigma(1, 0, R_gold, R_gold), rho_object) # mean-free path length of object (nm)
l_l = l(A_water, sigma(2, 1, R_H, R_O), rho_water) # mean-free path length of liquid (nm)
l_SiN = l(A_SiN, sigma(3, 4, R_Si, R_N), rho_SiN) # mean-free path length of SiN membrane (nm)

#print(l_o)
#print(l_l)
#print(l_SiN)

#print(((l_o*l_l)/(l_o-l_l)))

X = ((l_o*l_l)/(l_o-l_l))
Y = (1/l_l)
Z = 2*((t_SiN)/(l_SiN))

print((-(H*Y)-Z))
print(N_0)
print(((-(H*Y)-Z)/(N_0)))
print((np.exp(-(H*Y)-(Z))/(N_0)))

def d_SNR(H, t_SiN):
    #d_SNR = ((l_o*l_l)/(l_o-l_l))*(np.log(-3*(((np.exp(-(H/l_l)-((2*t_SiN)/l_SiN)))/N_0)**(1/2))+np.exp(-(H/l_l)-((2*t_SiN)/l_SiN)))+(H/l_l)+((2*t_SiN)/l_SiN))
    d_SNR = X*(np.log((-3*(np.exp(-(H*Y)-(Z))/N_0))+(np.exp(-(H*Y)-(Z))))+(H*Y)+(Z))
    return d_SNR
    
def d_blur(H, t_SiN):
    d_blur = 1.3*(((Lambda)**2)/(2*np.pi*a_h))*(H**(1.5))*(((rho_water*N_A)/(3*np.pi*A_gold))**(1/2))*Z_gold*(1 + (E/E_0))
    return d_blur
    
def d_chromatic_aberration(H, t_SiN):
    d_chromatic_aberration = C_c*alpha_0*((((N_A*(e**4)*rho_water*H)/(2*np.pi*(epsilon_0**2)*((v/c)**2))*(Z_gold/A_gold))/(2*E)))*((1+(E/E_0))/(1+(E/(2*E_0))))
    return d_chromatic_aberration

def d_OR(H, t_SiN):
    d_OR = ((d_SNR)**2 + (d_chromatic_aberration)**2 + (d_blur)**2)**(1/2)
    return d_OR 




#d_scherzer = (3/4)**(1/4)*(C_s)**(1/4)*(Lambda)**(3/4) 

d_OR_20 = d_OR(H, 20)
d_OR_40 = d_OR(H, 40)
d_OR_60 = d_OR(H, 60)
d_OR_80 = d_OR(H, 80)
d_OR_100 = d_OR(H, 100)

plt.figure(1)
plt.plot(H, d_OR_20, label = "10 nm Membrane")
plt.plot(H, d_OR_40, label = "20 nm Membrane")
plt.plot(H, d_OR_60, label = "30 nm Membrane")
plt.plot(H, d_OR_80, label = "40 nm Membrane")
plt.plot(H, d_OR_100, label = "50 nm Membrane")
plt.ylabel("d" + r"$\_{OR}$" + " (nm)")
#plt.xlim([0*10**(-7), 2000*10**(-7)]) #need to convert from cm to nm
#plt.ylim([0, 20])
plt.xlabel("Nominal Liquid Thickness (nm)")
plt.title("Resolution of LCTEM")
plt.legend()
