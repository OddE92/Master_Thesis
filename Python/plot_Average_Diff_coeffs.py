import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

filenameBase = "Data/SSH/"
filenameLM = ""
filenameE = ""
filenameEigenvalue= "eigenvalues.dat"

numRanks = 10
numParticlePerRank = 100

LM = 0
r_vect_len = 0
eigenval_vect = np.zeros(12)
i = 0

for filenameLM in ["LM10/", "LM150/"]:

    for j in range(3, 9):
        
        filenameE = "Ee1" + str(j) + "_"

        if j > 8: 
            print("filenameE out of range")
            break

        eigenval_current = np.loadtxt(filenameBase + filenameE + filenameLM + filenameEigenvalue) 
        #print(np.shape(eigenval_current))

        eigenval_len = len(eigenval_current)

        eigenval_vect[i * 6 + j-3] = eigenval_current[eigenval_len-1, 3]
        #print(filenameE + " eigenvalue: " + str(eigenval_vect[i * 6 + j-3]))

    i += 1

R_l = np.zeros(6)
for i in range(13, 19):
    R_l[i-13] = 0.270252 * 10**(i-15)
#print(R_l)

#region Calculate D_iso for LM = 10/150pc

g = 5.0/3.0
L_min = 0.027
L_max1 = 10.0
L_max2 = 150.0

L_maxmin1 = L_min/L_max1
L_maxmin2 = L_min/L_max2
L_c1 = 0.5 * L_max1 * (g-1)/g * ( (1-L_maxmin1)**g / (1-L_maxmin1)**(g-1)  )                #correlation length from Kristians paper
L_c2 = 0.5 * L_max2 * (g-1)/g * ( (1-L_maxmin2)**g / (1-L_maxmin2)**(g-1)  )                #Charged particle movement in turbulent MF

L01 = L_c1 / (2 * np.pi)                      #L_coh / 2*pi
L02 = L_c2 / (2 * np.pi) 
c = 0.3064                                    #pc/year
ccm = 2.998e10

res       = 50
E_iso     = np.logspace(12, 19, res)
R_l_iso   = np.logspace(np.log10(R_l[0]/10), np.log10(R_l[len(R_l)-1]*10), res)
D_iso_10  = np.zeros(res)
D_iso_150 = np.zeros(res)

unit_coeff = 3.019 * (10 ** 29)

D_iso_10[:]  = unit_coeff * c * (L01 / 3.0) *( (R_l_iso[:]/L01)**(2 - g) + (R_l_iso[:]/L01)**2 )
D_iso_150[:] = unit_coeff * c * (L02 / 3.0) *( (R_l_iso[:]/L02)**(2 - g) + (R_l_iso[:]/L02)**2 )

#endregion

E = [10**14, 10**15, 10**16, 5*10**16, 10**17, 5*10**17, 10**18, 5*10**18]
E0 = [10**13, 10**14, 10**15, 10**16, 10**17, 10**18]

#############################################################################################################################

figure = plt.figure(1)
ax = plt.gca()

#print len(E0)
#print len(R_l)

plt.plot(E0, unit_coeff * eigenval_vect[0:len(R_l)], '-s', color='red')
plt.plot(E0, unit_coeff * eigenval_vect[len(R_l)::],  '-o', color='green')

#print(E_iso)

plt.plot(E_iso, D_iso_10, '--', color='red')
plt.plot(E_iso, D_iso_150, '--', color='green')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title('Diffusion Tensor after transition time')
ax.set_ylabel('D $[\mathrm{cm}^2/\mathrm{s}]$')
ax.set_xlabel('Energy [eV]')
plt.legend(('$L_{max} = 10$pc', '$L_{max} = 150$pc', '$D_\mathrm{iso}(L_{max} = 10\mathrm{pc})$', '$D_\mathrm{iso}(L_{max} = 150\mathrm{pc})$'))

plt.show()


###############################################################################################################################

#region old code

# #D0 = np.loadtxt('data/kcor_eigenvalues_Ee14_LM150_nk100_N100.dat')
# D1 = np.loadtxt('data/Ee15_LM150/eigenvalues.dat')
# D2 = np.loadtxt('data/Ee16_LM150/eigenvalues.dat') 
# #D3 = np.loadtxt('data/kcor_eigenvalues_E5e16_LM150_nk500_N100.dat') 
# D4 = np.loadtxt('data/Ee17_LM150/eigenvalues.dat') 
# #D5 = np.loadtxt('data/kcor_eigenvalues_E5e17_LM150_nk500_N100.dat') 
# D6 = np.loadtxt('data/Ee18_LM150/eigenvalues.dat') 
# #D7 = np.loadtxt('data/kcor_eigenvalues_E5e18_LM150_nk500_N100.dat') 

# #D15 = np.loadtxt('data/kcor_eigenvalues_Ee14_LM10_nk100_N100.dat')
# D8  = np.loadtxt('data/Ee15_LM10/eigenvalues.dat') 
# D9  = np.loadtxt('data/Ee16_LM10/eigenvalues.dat')
# #D10 = np.loadtxt('data/kcor_eigenvalues_E5e16_LM10_nk500_N100.dat') 
# D11 = np.loadtxt('data/Ee17_LM10/eigenvalues.dat') 
# #D12 = np.loadtxt('data/kcor_eigenvalues_E5e17_LM10_nk500_N100.dat') 
# D13 = np.loadtxt('data/Ee18_LM10/eigenvalues.dat') 
# #D14 = np.loadtxt('data/kcor_eigenvalues_E5e18_LM10_nk500_N100.dat') 

# R_l0 = 0.0270252
# R_l1 = 0.270252                             #Calculate with calculate_trajectory.exe
# R_l2 = 10*R_l1
# R_l3 = 10*R_l2
# R_l4 = 10*R_l3

# #D_iso0 = 3.019 * (10 ** 29) * (c * L01 / 3.0)*( (R_l0 / L01)**(2 - 5/3) + (R_l0 / L01)**2 )      # cm^2/s
# D_iso1 = 3.019 * (10 ** 29) * (c * L01 / 3.0)*( (R_l1 / L01)**(2 - 5/3) + (R_l1 / L01)**2 )      # cm^2/s
# D_iso2 = 3.019 * (10 ** 29) * (c * L01 / 3.0)*( (R_l2 / L01)**(2 - 5/3) + (R_l2 / L01)**2 )      # cm^2/s
# D_iso3 = 3.019 * (10 ** 29) * (c * L01 / 3.0)*( (R_l3 / L01)**(2 - 5/3) + (R_l3 / L01)**2 )      # cm^2/s
# D_iso4 = 3.019 * (10 ** 29) * (c * L01 / 3.0)*( (R_l4 / L01)**(2 - 5/3) + (R_l4 / L01)**2 )      # cm^2/s

# #D_iso00 = 3.019 * (10 ** 29) * (c * L02 / 3.0)*( (R_l0 / L02)**(2 - 5/3) + (R_l0 / L02)**2 )      # cm^2/s
# D_iso10 = 3.019 * (10 ** 29) * (c * L02 / 3.0)*( (R_l1 / L02)**(2 - 5/3) + (R_l1 / L02)**2 )      # cm^2/s
# D_iso20 = 3.019 * (10 ** 29) * (c * L02 / 3.0)*( (R_l2 / L02)**(2 - 5/3) + (R_l2 / L02)**2 )      # cm^2/s
# D_iso30 = 3.019 * (10 ** 29) * (c * L02 / 3.0)*( (R_l3 / L02)**(2 - 5/3) + (R_l3 / L02)**2 )      # cm^2/s
# D_iso40 = 3.019 * (10 ** 29) * (c * L02 / 3.0)*( (R_l4 / L02)**(2 - 5/3) + (R_l4 / L02)**2 )      # cm^2/s

# D150_plot = [scal*D1[45, 3], scal*D2[45, 3], scal*D4[54, 3], scal*D6[63, 3]]
# D10_plot  = [scal*D8[45, 3], scal*D9[45, 3], scal*D11[54, 3], scal*D13[63, 3]]


#endregion