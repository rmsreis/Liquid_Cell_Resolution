import numpy as np
import matplotlib.pyplot as plt

c = 3*(10**8) # speed of light (m/s)
E_0 = 511*(10**3) # electron rest energy (eV)
N_A = 6.022*(10**23) # Avagadro's number
e = 1.602*(10**(-19)) # charge of electron (C)
a_h = 5.29*(10**(-11)) # bohr radius (m)
#epsilon_0 = 8.854*10**(-21) # permittivity of free space (C^2 g^-1 cm^-3 s^2)
epsilon_0 = 8.854*10**(-12) # permittivity of free space (C^2 N^-1 m^-2)

E = 300*(10**3) # incident electron energy (eV)
Lambda = 1.9687*(10**(-12)) # relativistic wavelength of electron (m) https://www.jeol.com/words/emterms/20121023.071258.php#gsc.tab=0
C_c = 2.4*(10**(-3)) # chromatic aberration coefficient (m)
C_s = 1.6*(10**(-3)) # spherical aberration coefficient (m)
d_objective_aperture = 60*(10**(-6)) # diameter of objective aperture (m)


alpha_0 = np.arcsin(((Lambda)/(d_objective_aperture))) # half-angle of objective aperture (RETURNS ZERO)

print((Lambda/(d_objective_aperture)))
print(alpha_0)

v = c*((1-(1/((1+(E/E_0))**2)))**(1/2)) # relativistic velocity (m/s)

H = np.arange(0*10**(-7), 2100*10**(-7), 100*10**(-7)) #liquid thickness (cm) #deltaE should have H in cm, but why?

t_SiN = 20*10**(-9) # membrane thickness (m)


flux = 0.1 # dose rate (electrons/Angstrom^2/s)
exposure_time = 1 # exposure time (s)
probe_diameter = 12.59*10**4 # beam/probe diameter (Angstrom)
probe_radius = probe_diameter/2 # beam/probe radius (Angstrom)
probe_area = np.pi*((probe_radius)**2) # beam/probe area (Angstrom^2)
N_0 = flux*exposure_time*probe_area # number of incident electrons 

deltaEcoeff = ((N_A*(e**4))/(2*np.pi*(epsilon_0**2)))
print("deltaE coefficient = " + str(deltaEcoeff))

def sigma(x, y, sigma_X, sigma_Y):
    sigma = ((x/(x+y))*sigma_X) + ((y/(x+y))*sigma_Y)
    return sigma

def sigma2(Z):
    numerator = (Z**2)*((a_h*(Z**(-1/3)))**2)*(Lambda**2)*((1+(E/E_0))**2)
    beta_omega = ((2*np.pi*a_h*alpha_0)/(Lambda*(Z**(1/3))))
    denominator = (np.pi*(a_h**2))*(1+(beta_omega**2))
    sigma2 = ((numerator)/(denominator))
    return sigma2

def l(A, sigma, rho):
    l = A/(sigma*rho*N_A)
    return l


A_water = ((2*1.008)+(1*15.999)) # atomic weight of water (amu)
A_SiN = ((3*28.084)+(4*14.0064)) # atomic weight of SiN (amu)
rho_water = 1.0*10**6 # specimen density of water (g/m^3)
#rho_water = 1.0 # specimen density of water (g/cm^3)
rho_SiN = 3.17*10**6 # specimen density of SiN (g/m^3)
#rho_SiN = 3.17 # specimen density of SiN (g/cm^3)

Z_H = 1 # atomic number of Hydrogen
Z_C = 6 # atomic number of Carbon
Z_N = 7 # atomic number of Nitrogen
Z_O = 8 # atomic number of Oxygen
Z_Si = 14 # atomic number of Silicon

Z_gold = 79 # atomic number of Gold
A_gold = 196.967 # atomic weight of Gold (amu)
rho_object = 19.3*10**6 # specimen density of Gold (g/m^3)
#rho_object = 19.3 # specimen density of Gold (g/cm^3)
Z_water = 4.7 # average atomic number for water ???

l_o = l(A_gold, sigma(1, 1, sigma2(Z_gold), sigma2(Z_gold)), rho_object) # mean-free path length of object (nm)
l_l = l(A_water, sigma(2, 1, sigma2(Z_H), sigma2(Z_O)), rho_water) # mean-free path length of liquid (nm)
l_SiN = l(A_SiN, sigma(3, 4, sigma2(Z_Si), sigma2(Z_N)), rho_SiN) # mean-free path length of SiN membrane (nm)

print(str(sigma(1, 1, sigma2(Z_gold), sigma2(Z_gold))))
print(str(sigma(2, 1, sigma2(Z_H), sigma2(Z_O))))
print(str(sigma(3, 4, sigma2(Z_Si), sigma2(Z_N))))
print("l_o = " + str(l_o))
print("l_l = " + str(l_l))
print("l_SiN = " + str(l_SiN))

#print(((l_o*l_l)/(l_o-l_l)))

X = ((l_o*l_l)/(l_o-l_l))
Y = (1/l_l)
W = (2/(l_SiN))

print((-(H*Y)-(t_SiN*W)))
print(N_0)
print(((-(H*Y)-(t_SiN*W))/(N_0)))
print((np.exp(-(H*Y)-(t_SiN*W))/(N_0)))

def d_SNR(H, t_SiN):
    #d_SNR = ((l_o*l_l)/(l_o-l_l))*(np.log(-3*(((np.exp(-(H/l_l)-((2*t_SiN)/l_SiN)))/N_0)**(1/2))+np.exp(-(H/l_l)-((2*t_SiN)/l_SiN)))+(H/l_l)+((2*t_SiN)/l_SiN))
    d_SNR = X*(np.log((-3*(np.exp(-(H*Y)-(t_SiN*W))/N_0))+(np.exp(-(H*Y)-(t_SiN*W))))+(H*Y)+(t_SiN*W))
    return d_SNR
    
def d_blur(H, t_SiN):
    d_blur = 1.3*(((Lambda)**2)/(2*np.pi*a_h))*(H**(1.5))*(((rho_water*N_A)/(3*np.pi*A_water))**(1/2))*(Z_water*(1 + (E/E_0)))
    return d_blur
    
def d_chromatic_aberration(H, t_SiN):
    d_chromatic_aberration = C_c*alpha_0*((((N_A*(e**4)*rho_water*H)/(2*np.pi*(epsilon_0**2)*((v/c)**2))*(Z_water/A_water))/(2*E)))*((1+(E/E_0))/(1+(E/(2*E_0))))
    return d_chromatic_aberration

def d_OR(H, t_SiN):
    d_OR = ((d_SNR(H, t_SiN))**2 + (d_chromatic_aberration(H, t_SiN))**2 + (d_blur(H, t_SiN))**2)**(1/2)
    return d_OR 




#d_scherzer = (3/4)**(1/4)*(C_s)**(1/4)*(Lambda)**(3/4) 

d_SNR_20 = d_SNR(H, 20)
d_SNR_40 = d_SNR(H, 40)
d_SNR_60 = d_SNR(H, 60)
d_SNR_80 = d_SNR(H, 80)
d_SNR_100 = d_SNR(H, 100)

d_blur_20 = d_blur(H, 20)
d_blur_40 = d_blur(H, 40)
d_blur_60 = d_blur(H, 60)
d_blur_80 = d_blur(H, 80)
d_blur_100 = d_blur(H, 100)

d_chromatic_aberration_20 = d_chromatic_aberration(H, 20)
d_chromatic_aberration_40 = d_chromatic_aberration(H, 40)
d_chromatic_aberration_60 = d_chromatic_aberration(H, 60)
d_chromatic_aberration_80 = d_chromatic_aberration(H, 80)
d_chromatic_aberration_100 = d_chromatic_aberration(H, 100)

d_OR_20 = d_OR(H, 20)
d_OR_40 = d_OR(H, 40)
d_OR_60 = d_OR(H, 60)
d_OR_80 = d_OR(H, 80)
d_OR_100 = d_OR(H, 100)


plt.figure(1)
plt.plot(H, d_SNR_20, label = "10 nm Membrane")
plt.plot(H, d_SNR_40, label = "20 nm Membrane")
plt.plot(H, d_SNR_60, label = "30 nm Membrane")
plt.plot(H, d_SNR_80, label = "40 nm Membrane")
plt.plot(H, d_SNR_100, label = "50 nm Membrane")
plt.ylabel("d" + r"$\d_{SNR}$" + " (nm)")
#plt.xlim([0*10**(-7), 2000*10**(-7)]) #need to convert from cm to nm
#plt.ylim([0, 20])
#plt.xticks(np.arange(0, 2100, 100))
plt.xlabel("Nominal Liquid Thickness (nm)")
plt.title("Resolution of LCTEM")
plt.legend()

plt.figure(2)
plt.plot(H, d_blur_20, label = "10 nm Membrane")
plt.plot(H, d_blur_40, label = "20 nm Membrane")
plt.plot(H, d_blur_60, label = "30 nm Membrane")
plt.plot(H, d_blur_80, label = "40 nm Membrane")
plt.plot(H, d_blur_100, label = "50 nm Membrane")
plt.ylabel("d" + r"$\d_{blur}$" + " (nm)")
#plt.xlim([0*10**(-7), 2000*10**(-7)]) #need to convert from cm to nm
#plt.ylim([0, 20])
#plt.xticks(np.arange(0, 2100, 100))
plt.xlabel("Nominal Liquid Thickness (nm)")
plt.title("Resolution of LCTEM")
plt.legend()

plt.figure(3)
plt.plot(H, d_chromatic_aberration_20, label = "10 nm Membrane")
plt.plot(H, d_chromatic_aberration_40, label = "20 nm Membrane")
plt.plot(H, d_chromatic_aberration_60, label = "30 nm Membrane")
plt.plot(H, d_chromatic_aberration_80, label = "40 nm Membrane")
plt.plot(H, d_chromatic_aberration_100, label = "50 nm Membrane")
plt.ylabel("d" + r"$\d_{chromatic aberration}$" + " (nm)")
#plt.xticks(np.arange(0, 2100, 100))
#plt.xlim([0*10**(-7), 2000*10**(-7)]) #need to convert from cm to nm
#plt.ylim([0, 20])
plt.xlabel("Nominal Liquid Thickness (nm)")
plt.title("Resolution of LCTEM")
plt.legend()



plt.figure(2)
plt.plot(H, d_OR_20, label = "10 nm Membrane")
plt.plot(H, d_OR_40, label = "20 nm Membrane")
plt.plot(H, d_OR_60, label = "30 nm Membrane")
plt.plot(H, d_OR_80, label = "40 nm Membrane")
plt.plot(H, d_OR_100, label = "50 nm Membrane")
plt.ylabel("d" + r"$\d_{OR}$" + " (nm)")
#plt.xticks(np.arange(0, 2100, 100))
#plt.xlim([0*10**(-7), 2000*10**(-7)]) #need to convert from cm to nm
#plt.ylim([0, 20])
plt.xlabel("Nominal Liquid Thickness (nm)")
plt.title("Resolution of LCTEM")
plt.legend()
