import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

def HTC_boil(dT): #heat transfer coefficient for boiling in function of temperature difference dT, estimated from chart, W/m2K
    htc = 152384.89074111 * dT ** 3 - 29229.8943538828 * dT ** 2 + 32570.6912303901 * dT - 690.763698082839
    return htc

def conduction_SS(T): #heat conduction coefficient of stainless steel in function of temperature T, estimated from chart, W/mK
    lam = -0.000000000000000486388582907927*T**6 + 0.000000000001885468436921880*T**5 - 0.00000000282179007851548*T**4 + 0.00000205319853863096*T**3 - 0.000758862385361664*T**2 + 0.152496009044047*T - 0.253795584108731
    return lam

def conduction_Cu(T): #heat conduction coefficient of copper in function of temperature T, estimated from chart, W/mK
    lam = -0.0117624718735243*T**3 + 0.0925982341332305*T**2 + 41.4673043235424*T + 1.35142280296072
    return lam

#Liquid helium flows through pipes and forced convection occurs. In case of helium, refprop should be used because the fact that coolprop isn't sufficiently precise:
def HTC_forced_conv_He(D, m, T): #(tube diameter, mass flow, temperature) heat transfer coefficient (forced convection) of helium, W/m2K
    Pr = 0.270854743014979*T**3 - 3.33975363875197*T**2 + 13.9239393304794*T - 18.9348695818743 #Prandtl number from refprop, -
    lam = 0.00021797032525761*T**3 - 0.00392564761753178*T**2 + 0.0224233218358592*T - 0.0210118147354211 #heat conductivity from refprop, W/mK
    kin_visc = 0.00000457636892847724*T**3 - 0.0000544948324032735*T**2 + 0.000203235396034521*T + 0.0000388609304085505 #kinematic viscosity from refprop, m2/s
    ro= -4.55149134862586*T**3 + 53.2457533093111*T**2 - 221.498152095213*T + 462.601030218146 #density, kg/m3
    U = m*ro/(3.14/4*D**2) #velocity in the pipe, m/s
    Re = U*D/kin_visc #Reynolds number of helium, -
    #Generally, Nusselt number Nu depends on Reynolds number:
    if Re <=10000 and Re > 0:
        Nu = 4.36
    elif Re > 10000:
        Nu = 0.023*Re**0.8*Pr*0.3 #correlation for Nu from Cengel heat transfer approach
    else:
        print("Reynolds number below zero!")
    htc = Nu*lam/D #heat transfer coefficient (forced convection) of helium, W/m2K
    #print(Nu)
    #print(lam)
    #print(Re)
    return htc

#Pressure drops resulting from friction during fluid flow:
def dp_He(D, m, T): #(tube diameter, mass flow, temperature)
    kin_visc = 0.00000457636892847724 * T ** 3 - 0.0000544948324032735 * T ** 2 + 0.000203235396034521 * T + 0.0000388609304085505  #kinematic viscosity from refprop, m2/s
    ro = -4.55149134862586 * T ** 3 + 53.2457533093111 * T ** 2 - 221.498152095213 * T + 462.601030218146  #density, kg/m3
    U = m / ro / (4 / 3.14 * D ** 2)  #velocity in the pipe, m/s
    Re = U * D / kin_visc  #Reynolds number, -
    f = (0.79*np.log(Re)-1.64)**(-2) #friction factor - previously flow character was checked (turbulent flow) so I will use formula for turbulent flow taken from Cengel
    dp = f*1./D*ro*U**2/2 #presser drops - 1/D instead of L/D because we want to calculate pressure drops per 1m of a pipe. Pressure drop on whole pipe will be calculated later
    return dp

def Main(D_tube, m_He, material, T_He_bath, T_He_coil_in,p_He_coil_in, g_tube, Q_duty): #(tube diameter, mass flow, material, helium bath temperature, helium temperature at inlet, helium pressure at inlet, pipe thickness, heat duty)
    h_He_coil_in = CP.PropsSI('H','T',T_He_coil_in,'P',p_He_coil_in,'Helium')       #enthalpy at inlet
    h_He_coil_out = h_He_coil_in - Q_duty/m_He                                      #enthalpy at outlet
    T_He_coil_out = CP.PropsSI('T','H',h_He_coil_out,'P',p_He_coil_in,'Helium')     #temperature at outlet, K
    MTD = ((T_He_coil_in-T_He_bath)+(T_He_coil_out-T_He_bath))/2
    T_He_coil = (T_He_coil_in +T_He_coil_out)/2
    dT = 0.025
    Convergence = False
    T_wall = (T_He_bath + T_He_coil)/2
    i = 0
    while Convergence == False:
        i = i +1
        #print(i)
        #print(dT)
        htc_fc = HTC_forced_conv_He(D_tube, m_He, T_He_coil )
        htc_boil = HTC_boil(dT)
        htc_wall = 0
        if material == 'Cu':
            htc_wall = conduction_Cu(T_wall)/g_tube
        elif material == 'ss':
            htc_wall = conduction_SS(T_wall)/g_tube
        else:
            print("There is no such material specified. Choose between 'Cu' or 'ss'")
        ohtc = 1/(1/htc_fc + 1/htc_boil + 1/htc_wall)               #overall heat transfer coefficient, W/m2K
        q = ohtc * MTD                                              #MTD - mean temperature diff - we don't need to calculate logarithmic mean temerature difference because temperature drop is almost linear)
        dT_1 = q/htc_boil
        #print(dT_1)
        dT_fc = q/htc_fc                                            #forced convection temperature difference in order to determine wall temperature more precisely
        T_wall = ((T_He_bath + dT_1)+(T_He_coil - dT_fc))/2
        epsilon = abs(dT - dT_1)/dT_1
        #print(epsilon)
        if epsilon < 0.05:
            Convergence = True
        else:
            dT = dT + 0.0001
    #print(str(round(htc_fc, 2)) + " - heat transfer coefficient (forced convection), W/m2K")
    #print(str(round(htc_wall,2)) + " - heat transfer coefficient of wall, W/m2K")
    #print(str(round(htc_boil,2)) + " - heat transfer coefficient for boiling, W/m2K")
    #print(str(round(ohtc,2)) + " - overall heat transfer coefficient, W/m2K")
    A = Q_duty/q                                                    #required heat transfer area, m2
    L = A/(3.14 * (D_tube-2*g_tube))*1.2                            #where 1.2 is a safety factor
    dp_m = dp_He(D_tube, m_He, T_He_coil_in)
    dp = dp_m*L
    return [L, dp]

#print(Main(0.01,0.4, 'ss',4.6, 4.8, 300000, 0.001, 150)) #(D_tube, m_He, material, T_He_bath, T_He_coil_in,p_He_coil_in, g_tube, Q_duty)

D_list = []
L_list = []
dp_list = []
for i in range(4, 101):
    D_tube = i/1000.
    results = Main(D_tube,0.4, 'Cu',4.6, 4.8, 300000, 0.001, 150)
    L = results[0]
    dp = results[1]
    D_list.append(D_tube)
    L_list.append(L)
    dp_list.append(dp)
#print(D_list)
#print(L_list)

plt.figure(1)
plt.plot(D_list, L_list)
plt.xlabel('Tube diameter, m')
plt.ylabel('Length of the coil, m')

plt.figure(2)
plt.plot(D_list, dp_list)
plt.yscale(value = 'log')
plt.xlabel('Tube diameter, m')
plt.ylabel('Pressure drop, Pa')

plt.show()