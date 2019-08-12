# -*- coding: utf-8 -*-

#
# rocket-launch.py
# perform a launch for the water rocket
#
# hvf@hvf.ch 31.05.17 04.06.2017 06.06.2017

import pint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = (16.0, 9.0)

ureg = pint.UnitRegistry(system='mks')
ureg.setup_matplotlib()

#%% Constants

dair = 1.225 * ureg.kilogram / ureg.meter**3   # [kg/m^3] # air density at sea level, ISA

dwater = 1000. * ureg.kilogram / ureg.meter**3  # [kg/m^3] # water density

pAtm = 101325 * ureg.pascal  # atmospheric pressure

g = 9.80665 * ureg.meter / ureg.second**2     # earth gravitation [m/s2]
   
Rgas = 287.058 * ureg.joule / ureg.kilogram / ureg.kelvin # [J/kg/K] ideal gas constant divided by gas molecular weight

Temp = 300 * ureg.kelvin # [K] constant gas temperature (absolute scale: K or R)

gamma = 1.4 # adiabatic index of water


#%% psi to pascal convenience
def psi_to_pascal(psi):
    pa = 6894.7572931783 * psi
    return pa

#%% 
class bottle_rocket(object):
 #_____________________________________________________________________________   
    def __init__(self, p0, V_water, V_bottle, d_bottle, d_nozzle, m_bottle, Cdrag, dt = 0.001 * ureg.second):
        
        # Model parameters

        # time step
        self.dt = dt.to_base_units() # sec
        
        # initial pressure in air chamber (empty volume)
        self.p0 = p0.to_base_units()    # Pa

        
        # total volume of water rocket
        self.Vtot = V_bottle.to_base_units() # 1.5 liter
        
        self.Vwater = V_water.to_base_units() # 0.5 # l
        
        # empty (air) volume of water rocket; rest is filled with water
#        V0 = 0.5e-3     # empty volume, 0.5 liter
#        V0 = 0.4e-3     # empty volume, 0.4 liter
#        V0 = 0.5e-3     # empty volume, 0.5 liter
#        V0 = 0.6e-3     # empty volume, 0.6 liter
#        V0 = 0.75e-3     # empty volume, 0.75 liter
#        V0 = 0.825e-3     # empty volume, 0.825 liter
#        V0 = 0.90e-3     # empty volume, 0.9 liter
#        V0 = 1.00e-3     # empty volume, 1.0 liter
#        V0 = 0.5e-3     # empty volume, 0.5 liter
        
        self.Vtot = V_bottle.to_base_units() #* 1e-3 # litres to m^3
        self.Vwater = V_water.to_base_units() #* 1e-3 # litres to m^3
        self.V0 = self.Vtot - self.Vwater
        
        
        self.p0 = pAtm + p0
        
        # cross section of water rocket (cylinder of 10.5 cm diameter)
        self.A  = (d_bottle**2*np.pi/4).to_base_units()    # m2, diam 10.5 cm
        
        # cross section of bottle opening, 2.5 cm diameter
        self.Q  = (d_nozzle**2*np.pi/4).to_base_units()    # m2, diam 2.5 cm
        
        # arbitrary throttling factor such that water does not exit too fast
        # we aim at a propulsion time of about 0.3 sec for 1 liter of water
        #Q = Q / 1.7     # throttling, arbitrary throttling factor
        
        # mass of empty bottle
        self.M  = m_bottle.to_base_units() #0.040      # kg, 40 gr empty weight of rocket (bottle)
        
        # drag coefficient of bottle 
        self.Cdrag = Cdrag #0.6
        
        # some variables and initial values
        a0 = 0. * ureg.meter/ureg.second**2
        v0 = 0. * ureg.meter/ureg.second
        s = 0. * ureg.meter       # [m], height of rocket above ground
        text = ""
        
        # plotting arrays
        self.timelist = np.array([])
        self.heights =  np.array([])
        self.velocs =  np.array([])
        self.accels =  np.array([])
        self.jetspeed =  np.array([])
        
        
#        print ()
#        print ("  Vtot    V0      A        Q         weight  dwater   dair   Cdrag")
#        print ("  [m³]    [m³]    [m2]     [m2]      [kg]    [kg/m³]  [kg/m³]      ")
#        
#        print ("{7:7.2e} {0:7.2e} {1:7.2e} {2:7.2e} {3:7.3f}   {4:5.1f}  {5:7.4f} {6:4.2f}".format(self.V0.magnitude, self.A.magnitude, self.Q.magnitude, self.M.magnitude, dwater.magnitude, dair.magnitude, self.Cdrag, self.Vtot.magnitude))
#        
#        print ()
#        print (" t       p         V         v       a        ds    vj       dt  Thrust mdot  m   Drag brake")
#        print (" [s]     [Pa]      [m³]      [m/s]   [m/s²]   [m]   [m/s]    [s]  [N]   [kg]  [kg]  [N]   [m/s²]")
        
        p = self.p0
        V = self.V0
        a = a0
        v = v0
        t = 0.0 * ureg.second
        vj = 0. * ureg.meter/ureg.second
        
        #self.doPrintStep (t, p, V, v, a, s, vj, dt, text)
              
#        self.timelist = np.append(self.timelist, t)
#        self.heights = np.append(self.heights, s)
#        self.velocs = np.append(self.velocs, v0)
#        self.accels = np.append(self.accels, a0)
#        self.jetspeed = np.append(self.jetspeed, 0. * ureg.meter/ureg.second)
        
        # find out maximal values and keep
        self.hmax = 0. * ureg.meter; self.thmax = 0. * ureg.second
        self.vmax = 0. * ureg.meter/ureg.second; self.tvmax = 0. * ureg.second
        self.amax = 0. * ureg.meter/ureg.second**2; self.tamax = 0. * ureg.second
        self.jmax = 0. * ureg.meter/ureg.second; self.tjmax = 0. * ureg.second
        
        # let the rocket fly - max 1000 time steps 
        #for i in range(1000):
        while s >= 0. * ureg.meter:
            p, V, v, a, vj, ds, dt1, text = self.doTimeStep (dt, p, V, self.Vtot, v, a, self.A, self.Q, self.M, dwater, dair, self.Cdrag)
            t = t + dt1
            s = s + ds
            #self.doPrintStep (t, p, V, v, a, s, vj, dt1, text)
            
            self.timelist = np.append(self.timelist, t)
            self.heights = np.append(self.heights, s)
            self.velocs = np.append(self.velocs, v)
            self.accels = np.append(self.accels, a)
            self.jetspeed = np.append(self.jetspeed, vj)
            
            # find maximal values and times
            if (s > self.hmax): 
                self.hmax = s 
                self.thmax = t
            if (v > self.vmax): 
                self.vmax = v 
                self.tvmax = t
            if (a > self.amax): 
                self.amax = a 
                self.tamax = t
            if (vj > self.jmax): 
                self.jmax = vj 
                self.tjmax = t
                
            # end of thrust phase: larger time steps
            if (dt1 < self.dt or vj == 0.0 * ureg.meter/ureg.second):
                self.dt = 0.01 * ureg.second   #0.01 sec
                
            # stop if falling below floor
            if (s < 0.0 * ureg.meter): 
                break
                
            # increase time steps if falling
            if (v <= 0.0 * ureg.meter/ureg.second): 
                dt = 0.1 * ureg.second 
                
            # break after 3 secs
            #if (t >= 3.0): break
        
#        print ("Maximum altitude: {0:.1f} m at {1:.3f} s".format(self.hmax.magnitude, self.thmax.magnitude))
#        print ("Maximum velocity: {0:.1f} m/s at {1:.3f} s".format(self.vmax.magnitude, self.tvmax.magnitude))
#        print ("Maximum acceleration: {0:.1f} m/s² at {1:.3f} s".format(self.amax.magnitude, self.tamax.magnitude))
#        print ("Maximum jet speed: {0:.1f} m/s at {1:.3f} s".format(self.jmax.magnitude, self.tjmax.magnitude))
 
#_____________________________________________________________________________       
# doTimeStep - performs one time step of moving rocket
    def doTimeStep (self, dt, p0, V0, Vtot, v0, a0, A, Q, M, dwater, dair, Cdrag):
        """ calculate one time step for a water rocket.
        Input parameters:
        dt  time step [s]
        p0  initial gas pressure [Pa]
        V0  initial gas volume [m³]
        Vtot    total volume of rocket [m³] - Vtot-V0 is volume of water/propellant
        v0  initial vertical velocity of rocket [m/s]
        a0  initial acceleration of rocket [m/s2]
        A   cross section of rocket (skin thickness is zero) [m2]
        Q   nozzle cross section [m2]
        M   empty mass of rocket [kg]
        dwater  density of water/propellant [kg/m³]
        dair    density of ambient air [kg/m³]
        Cdrag   drag coefficient of rocket [-]
    
        Return values:
        p1  final gas pressure [Pa]
        V1  final gas volume [m³]
        v1  final vertical velocity of rocket [m/s]
        a1  final acceleration of rocket [m/s]
        vj  velocity of water exhaust jet [m/s]
        ds1 increase in vertical position [m]
        dt1 actually used time step 
        text    additional values in a printable string
    
        Intermediate values:
        m   mass of water/propellant, dwater*(Vtot-V0) [kg]
        T   thrust of water jet [N]
        Vdot    volume flow [m³/s]
        mdot    mass flow [kg/s]
        D   air drag [N]
        """    
        
        # the interior pressure, p0, must be bigger than pAtm
        # or no water is propelled out of the rocket
    
        text = ""
    
        # Drag of flying rocket
        D = - Cdrag / 2.0 * A * dair * v0 * np.abs(v0)
    
        # check if pressure has equalised
        if p0 <= pAtm: 
            vj = 0.0  * ureg.meter/ureg.second  # no more thrust
            T = 0.0 * ureg.newton
            V1 = V0
            if (V0 >= Vtot): 
                V1 = Vtot 
                V0 = Vtot
                
            Vdot = 0.0 * ureg.meter**3/ureg.second
            mdot = 0.0 * ureg.kilogram/ureg.second
            
            # mass of water/propellant - there may still be water but no more p
            m = (Vtot - V1) * dwater
            dm = 0.0 * ureg.kilogram
            p1 = pAtm
                
            a1 = (-g + D / (M + m)).to_base_units() 
            adrag = (D / (M + m)).to_base_units()      # accel/decel due to drag
            text = "{0:3.0f} {1:.3f} {2:.3f} {3:.1f} {4:.2f}".format(0, 0, (M+m).magnitude, D.magnitude, adrag.magnitude)
            
        else:
            # check if all water is out, then thrust is only due air pressure
            if V0 >= Vtot:
                # calculate vj from Bernoulli, 1/2 *dwater*vj**2 = p0-pAtm
                if p0 > pAtm: 
                    p1 = p0 * np.exp(-Q/Vtot * np.sqrt(Rgas*Temp) * dt)
                    dt0 = np.log(pAtm/p0) / (-Q/Vtot * np.sqrt(Rgas*Temp))
                    if dt0 > dt:
                        vj = (np.sqrt(2.*(p0 - pAtm)/dair)).to_base_units() 
           
                    else:
                        vj = 0.0 * ureg.meter/ureg.second
                else:
                    vj = 0.0 * ureg.meter/ureg.second
                   
                # volume flow, Vdot = Q*vj
                Vdot = Q * vj
    
                # mass flow mdot = Q*vj*dair
                mdot = Vdot * dair
                
            else:
                # calculate vj from Bernoulli, 1/2 *dwater*vj**2 = p0-pAtm
                vj = (np.sqrt(2*(p0 - pAtm)/dwater)).to_base_units() 
    
                # volume flow, Vdot = Q*vj
                Vdot = Q * vj
    
                # mass flow mdot = Q*vj*dwater
                mdot = Vdot * dwater
                
                # check if there is still some water
                # how long does it take (at current conditions) to empty bottle?
                tmax = (Vtot - V0)/Vdot
                #if (tmax < dt): 
                 #   dt = tmax   # adjust dt down
            
    
            # thrust = mdot * vj
            T = mdot * vj
    
            dm = mdot * dt
            dV = Vdot * dt
    
            # new volume
            V1 = V0 + dV
            if V1 >= Vtot:
                V1 = Vtot
                p1 = p0 * np.exp(-Q/Vtot * np.sqrt(Rgas*Temp) * dt)
                
            else:
                # new pressure, ideal gas
                p1 = p0 * (V0/V1)**gamma
            
            # mass of water/propellant
            m = (Vtot-V1)*dwater  
                
            adrag = D/(M+m)     # accel/decel due to drag
            
            # thrust, water loss, total mass, drag
            text = "{0:3.0f} {1:.3f} {2:.3f} {3:.1f} {4:.2f}".format(T.magnitude, dm.magnitude, (M+m).magnitude, D.magnitude, adrag.magnitude)
        
            # acceleration from forces
            a1 = ((T - (M+m)*g + D)/(M+m)).to_base_units() 
    
        # intergrate a for v and v for s - may be improved later
#        print(v0)
#        print(a1)
#        print(dt)
        
        v1 = (v0 + a1*dt).to_base_units() 
        ds = (v0+v1)*dt/2.
    
        return (p1, V1, v1, a1, vj, ds, dt, text)
    
    def plot_results(self, fig = 0):
        ###########################################################
        # make a plot with subplots for therelevant stuff
        
        fig = plt.figure(fig)#, figsize=(10,10))
        #fig.subplots_adjust(bottom=0.10, left=0.15, top = 0.95, right=0.95)
        plt.subplot(221)
        
        plt.title("Altitude, maximum {0:.1f} m".format(self.hmax.magnitude))
        plt.plot(self.timelist, self.heights)
        plt.xlabel("Time [s]")
        plt.ylabel("Altitude [m]")
        plt.grid()
        
        plt.subplot(222)
        plt.title("Velocity, maximum {0:.1f} m/s".format(self.vmax.magnitude))
        plt.plot(self.timelist, self.velocs)
        plt.xlabel("Time [s]")
        plt.ylabel("Velocity [m/s]")
        plt.grid()
        
        plt.subplot(223)
        plt.title("Acceleration")
        plt.title("Acceleration, maximum {0:.1f} m/s²".format(self.amax.magnitude))
        plt.plot(self.timelist, self.accels)
        plt.xlabel("Time [s]")
        plt.ylabel("Acceleration [m/s²]")
        plt.grid()
        
        plt.subplot(224)
        plt.title("Jet Speed")
        plt.title("Jet Speed, maximum {0:.1f} m/s".format(self.jmax.magnitude))
        plt.plot(self.timelist, self.jetspeed)
        plt.xlabel("Time [s]")
        plt.ylabel("Jetspeed [m/s]")
        plt.grid()
        
        # create results file name with pressure p0 and volume V0
        resname = "rocket-results-{0:d}-{1:d}.pdf".format(int(self.p0.magnitude), int(self.V0.magnitude*1e6))
        plt.savefig(resname)
        print ("Created PDF file ", resname)
        figManager = plt.get_current_fig_manager()
        # figManager.window.showMaximized()
        fig.tight_layout()#rect=[0, 0.03, 1, 0.95])
        plt.show()
        plt.savefig('Default_Flight.png')
        
        
#_____________________________________________________________________________
# doPrintStep - print lots of stuff for each step
    def doPrintStep (self, t, p, V, v, a, ds, vj, dt, text):
        print ("{0:5.3f} {1:8.1f} {2:7.2e} {3:7.2f}   {4:7.2f} {5:7.2e} {6:7.3e} {7:.3f} {8:s}".format(t.magnitude, p.magnitude, V.magnitude, v.magnitude, a.magnitude, ds.magnitude, vj.magnitude, dt.magnitude, text))
            
#%% main
if __name__=="__main__":
    
    p0 = 60 * ureg.psi
    V_water = 0.66 * ureg.litre
    V_bottle = 2 * ureg.litre
    d_bottle = 10.0 * ureg.centimeter 
    d_nozzle = 10 * ureg.millimeter 
    m_bottle = 100 * ureg.grams
    Cdrag = 1.0 # good bottle rocket
#    Cdrag = 0.75 # typical bottle rocket
    dt = 0.001 * ureg.seconds
    
    figure = -1
    
#%%   
    print('Basic rocket')
    br_1 = bottle_rocket(p0 =  p0, 
                         V_water = V_water, 
                         V_bottle = V_bottle, 
                         d_bottle = d_bottle, 
                         d_nozzle = d_nozzle, 
                         m_bottle = m_bottle, 
                         Cdrag = Cdrag, 
                         dt = dt)
    
    figure+=1
    br_1.plot_results(figure)
    
#%% Optimise water volume
    
    print('Optimising water volume')
    maxs = np.array([])
    vols = np.array([])
    water_step = 0.025 * ureg.litre
    upper = (V_bottle + water_step).magnitude
    for v_water in np.arange(0, upper, water_step.magnitude):
        
        V_water = v_water * ureg.litre
        print(V_water)
        br_1 = bottle_rocket(p0 =  p0, 
                             V_water = V_water, 
                             V_bottle = V_bottle, 
                             d_bottle = d_bottle, 
                             d_nozzle = d_nozzle, 
                             m_bottle = m_bottle, 
                             Cdrag = Cdrag, 
                             dt = dt)
        
        maxs = np.append(maxs, br_1.hmax)
        vols = np.append(vols, V_water)
    
    figure+=1
    plt.figure(figure)
    fig, ax = plt.subplots()
    ax.set_title('Optimal Water Volume')
    ax.yaxis.set_units(ureg.meter)
    ax.xaxis.set_units(ureg.litre)
    ax.set_xlabel('Water Volume [ℓ]')
    ax.set_ylabel('Max. Height [m]')
    
    ax.plot(vols, maxs, 'tab:blue')
    ax.axvline(vols[np.argmax(maxs)], color='tab:green')
    ax.grid()
    fig.tight_layout()
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    plt.show()
    plt.savefig('Optimum_water.png')
    
    print('Max. height:', np.max(maxs))
    print('Optimum fill volume:', vols[np.argmax(maxs)])
    
    V_water = 0.275 * ureg.litre
 
#%% vary Pressure
    print('Varying pressure')
    figure+=1
    plt.figure(figure)
    fig, ax = plt.subplots()
    ax.set_title('Varied Pressure')
    
    brs = []
    step = 10 * ureg.psi
    upper = (60 * ureg.psi + step).magnitude
    for p in np.arange(step.magnitude, upper, step.magnitude):
        
        p0 = p * ureg.psi
        print(p0)
        br_1 = bottle_rocket(p0 =  p0, 
                             V_water = V_water, 
                             V_bottle = V_bottle, 
                             d_bottle = d_bottle, 
                             d_nozzle = d_nozzle, 
                             m_bottle = m_bottle, 
                             Cdrag = Cdrag, 
                             dt = dt)
    
        
        ax.yaxis.set_units(ureg.meter)
        ax.xaxis.set_units(ureg.seconds)
        ax.plot(br_1.timelist, br_1.heights, label= str(p) + ' psi')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Height [m]')
         
    ax.legend()
    ax.grid()
    fig.tight_layout()
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    plt.show()
    plt.savefig('Pressure.png')
    
    p0 = 60 * ureg.psi
        
    
#%% vary Water volume
    print('Varying water volume')
    figure+=1
    plt.figure(figure)
    fig, ax = plt.subplots()
    ax.set_title('Varied Water Volume')
    
    brs = []
    water_step = 0.25 * ureg.litre
    upper = (V_bottle / 2.0 + water_step).magnitude
    for v_water in np.arange(0, upper, water_step.magnitude):
        
        V_water = v_water * ureg.litre
        print(V_water)
        
        br_1 = bottle_rocket(p0 =  p0, 
                             V_water = V_water, 
                             V_bottle = V_bottle, 
                             d_bottle = d_bottle, 
                             d_nozzle = d_nozzle, 
                             m_bottle = m_bottle, 
                             Cdrag = Cdrag, 
                             dt = dt)
    
        
        ax.yaxis.set_units(ureg.meter)
        ax.xaxis.set_units(ureg.seconds)
        ax.plot(br_1.timelist, br_1.heights, label= str((br_1.Vwater.to(ureg.liter)).magnitude) + ' ℓ')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Height [m]')
         
    ax.legend()
    ax.grid()
    fig.tight_layout()
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    plt.show()
    plt.savefig('Water.png')
    
    V_water = 0.275 * ureg.litre
    
    
#%% vary Drag
    print('Varying drag coefficient')
    figure+=1
    plt.figure(figure)
    fig, ax = plt.subplots()
    ax.set_title('Varied Drag Coefficient')
    
    brs = []
    step = 0.2 
    upper = (1 + step)
    for Cdrag in np.arange(step, upper, step):
        print(Cdrag)
        br_1 = bottle_rocket(p0 =  p0, 
                             V_water = V_water, 
                             V_bottle = V_bottle, 
                             d_bottle = d_bottle, 
                             d_nozzle = d_nozzle, 
                             m_bottle = m_bottle, 
                             Cdrag = Cdrag, 
                             dt = dt)
    
        
        ax.yaxis.set_units(ureg.meter)
        ax.xaxis.set_units(ureg.seconds)
        ax.plot(br_1.timelist, br_1.heights, label= str(round(br_1.Cdrag, 1)))
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Height [m]')
         
    ax.legend()
    ax.grid()
    fig.tight_layout()
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    plt.show()
    plt.savefig('Drag.png')
    
    Cdrag = 0.35 # good bottle rocket
