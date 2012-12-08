#Library of functions called by SimpleBuildingEngine
import pandas as pd
import numpy as np

def WALLS(Btest=None):
    #Building height
    h_building = 2.7#[m]
    h_m_building = h_building / 2
    h_cl = 2.7# heigth of a storey

    #number of walls
    n_walls = 7
    A_fl = 48

    #WALLS CHARACTERISTICS
    #Orientation
    ori = pd.Series([('S'), ('W'), ('N'), ('E'), ('R'), ('F'), ('C')])
    #Surface azimuth
    surf_az = pd.Series([0, 90, 180 - 90, 0, 0, 0])
    #Slopes (90:vertical; 0:horizontal)
    slope = pd.Series([90, 90, 90, 90, 0, 0, 0])
    #Masks
    f_low_diff = pd.Series([1, 1, 1, 1, 1, 1, 1])
    f_low_dir = pd.Series([1, 1, 1, 1, 1, 1, 1])
    #U VALUES
    U_hopw = pd.Series([0.5144, 0.5144, 0.5144, 0.5144, 0.3177, 0, 0])
    U_lopw = pd.Series([3, 3, 3, 3, 3, 3, 3])
    U_fr = pd.Series([2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4])
    U_gl = pd.Series([3, 3, 3, 3, 3, 3, 3])


    if (Btest == 195 or Btest == 395):
        #SURFACES
        #Heavy Opaque walls
        A_hopw = pd.Series([21.6, 16.2, 21.6, 16.2, 48, 48, 48])
        #Windows
        A_wd = pd.Series([0, 0, 0, 0, 0, 0, 0])
        #Frame
        FWR = pd.Series([0, 0, 0, 0, 0, 0, 0])
        A_fr = FWR * A_wd
        #Glazing
        A_gl = A_wd - A_fr
        #Light Opaque walls
        A_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])

    elif (Btest == 200 or Btest == 210 or Btest == 230 or Btest == 240 or Btest == 250 or Btest == 400 or Btest == 410
          or Btest == 420 or Btest == 430 or Btest == 800):
        #Heavy Opaque walls
        A_hopw = pd.Series([9.6, 16.2, 21.6, 16.2, 48, 48, 48])
        #Windows
        A_wd = pd.Series([0, 0, 0, 0, 0, 0, 0])
        #Frame
        FWR = pd.Series([0, 0, 0, 0, 0, 0, 0])
        A_fr = FWR * A_wd
        #Glazing
        A_gl = A_wd - A_fr
        #Light Opaque walls
        A_lopw = pd.Series([12, 0, 0, 0, 0, 0, 0])

    elif (Btest == 270 or Btest == 320 or Btest == 600 or Btest == 640 or Btest == 650 or Btest == 810 or Btest == 900
          or Btest == 940 or Btest == 950 or Btest == 6001 or Btest == 9001 or Btest == 6501 or Btest == 9501):
        #Heavy Opaque walls
        A_hopw = pd.Series([9.6, 16.2, 21.6, 16.2, 48, 48, 48])
        #Windows
        A_wd = pd.Series([12, 0, 0, 0, 0, 0, 0])
        #Frame
        FWR = pd.Series([0, 0, 0, 0, 0, 0, 0])
        A_fr = FWR * A_wd
        #Glazing
        A_gl = A_wd - A_fr
        #Light Opaque walls
        A_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])

    elif (Btest == 300 or Btest == 620 or Btest == 920):
        #Heavy Opaque walls
        A_hopw = pd.Series([9.6, 16.2, 21.6, 16.2, 48, 48, 48])
        #Windows
        A_wd = pd.Series([0, 6, 0, 6, 0, 0, 0])
        #Frame
        FWR = pd.Series([0, 0, 0, 0, 0, 0, 0])
        A_fr = FWR * A_wd
        #Glazing
        A_gl = A_wd - A_fr
        #Light Opaque walls
        A_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])

    #Total    
    A_hopw_t = A_hopw.sum()
    A_wd_t = A_wd.sum()
    A_fr_t = A_fr.sum()
    A_lopw_t = A_lopw.sum()
    A_gl_t = max(0, A_wd_t - A_fr_t)
    A_t = A_hopw_t + A_lopw_t + A_wd_t + A_fr_t

    #CAPACITIES
    if (Btest == 800 or Btest == 900 or Btest == 920 or Btest == 940 or Btest == 950 or Btest == 9001 or Btest == 9501):
        C_hopw = ([145154, 145154, 145154, 145154, 18170, 112121, 0])
        C_lopw = ([0, 0, 0, 0, 0, 0, 0])
    else:
        C_hopw = ([14534, 14534, 14534, 14534, 18170, 19620, 0])
        C_lopw = ([0, 0, 0, 0, 0, 0, 0])
    
    C_m = sum((A_lopw * C_lopw + A_hopw * C_hopw))

    #Effective mass area [m^2]
    A_m = C_m ** 2 / sum((A_lopw * np.exp2(C_lopw) + A_hopw * np.exp2(C_hopw)))
    
    return n_walls, f_low_diff, f_low_dir, ori, surf_az, slope, A_t, A_fl, A_lopw_t, A_hopw_t, A_gl_t, A_fr_t, A_lopw,\
           A_hopw, A_gl, h_cl, C_m, A_m, U_hopw, U_lopw, U_fr, U_gl

def w_t_RH(p_atm=None, t=None, RH=None):
    from math import exp
    #Humidity ratio as function of drybulb temperature and humidity ratio
    p_w_s = exp((17.438 * t / (239.78 + t)) + 6.4147)#partial pressure of saturated water vapor
    p_w = RH * p_w_s
    w = (p_w * 0.62198) / (p_atm - p_w)
    return w

def ZENITHANG(Lat=None, Long=None, Long_st=None, n=None, h=None):
    from math import pi,cos,sin,acos
    from numpy import fix
    #ZENITH ANGLE
    #Ref: Duffie,J.A.,Beckman, W.A. 1980. Solar engineering of thermal
    #processes. 2nd Edition. John Wiley & Sons.
    #OUTPUTS
    # -h_sol: Solar time (in hours)
    # -h_sol_per: Solar time (in hours per day)
    # -phi: Latitude in radians
    # -delta: Declination angle in radians
    # -omega: Hour angle in radians
    # -theta_z: Zenith angle in radians, i.e. angle of incidence of beam radiation on a horizontal surface
    #INPUTS
    # -Lat: Latitude of the location (north positive) -90<Lat<90
    # -Long: Longitude of the location (west positive) 0<Long<180
    # -Long_st: Longitude of the standard meridian of the time zone
    # -n: day 1<n<365
    # -h: hour 1<h<8760

    #Angles in radians%
    phi = Lat * pi / 180

    #Summer time correction (Masy, 2008)
    epsilon_summer = 1

    #Equation of time (minutes)
    B = (n - 1) * 360 / 365 * pi / 180
    E = 229.2 * (0.000075 + 0.001868 * cos(B) - 0.032077 * sin(B) - 0.014615 * cos(2 * B) - 0.04089 * sin(2 * B))

    #Solar time (in hours)
    h_sol = h + (4 * (Long_st - Long) + E) / 60 - epsilon_summer

    #Solar time (in hours per day)
    h_sol_per_1 = h_sol - 24 * fix(h_sol / 24)
    if h_sol_per_1 <= 1E-6:
        h_sol_per = 24
    else:
        h_sol_per = h_sol_per_1

    #Declination (angular position of the sun at solar noon, north positive)
    #-23.45<delta<23.45
    delta = 23.45 * sin(360 * (284 + n) / 365 * pi / 180) * pi / 180#(daily basis, Cooper in Duffie & Beckmann)

    #Hour angle (morning negative, afternoon positive)
    omega = (h_sol_per - 12) * 15 * pi / 180

    #Zenith angle (between the vertical and the line to the sun)
    theta_z = max(1E-5, acos(cos(delta) * cos(phi) * cos(omega) + sin(delta) * sin(phi)))
    return phi, delta, omega, theta_z, h_sol

def CSITH(Lat=None, Long=None, Long_st=None, n=None, h=None):
    from math import cos,exp
    #Clear sky solar radiation
    #OUTPUTS
    # -I_th_cs: Clear sky theoretical solar radiation (in W/m2)
    #INPUTS
    # -Lat: Latitude of the location (north positive) -90<Lat<90
    # -Long: Longitude of the location (west positive) 0<Long<180
    # -Long_st: Longitude of the standard meridian of the time zone
    # -n: day 1<n<365
    # -h: hour 1<h<8760

    #Main angles and solar time for location
    phi, delta, omega, theta_z, h_sol = ZENITHANG(Lat, Long, Long_st, n, h)
    #Extraterrestrial radiation
    G_sc = 1353#W/m2 - Solar constant
    I_on = G_sc * (1 + 0.033 * cos(360 * (h_sol / 24) / 365))#Normal extraterrestrial radiation
    #Atmospheric transmittance for beam radiation (altitude = 0m)
    tau_b = 0.12814 + 0.7568875 * exp(-0.387225 / (cos(theta_z)))
    #Clear sky beam normal radiation
    I_cnb = I_on * tau_b
    #Clear sky horizontal beam radiation
    I_cb = I_cnb * cos(theta_z)
    #Atmospheric transmittance for diffuse radiation (altitude = 0m)
    tau_d = 0.271 - 0.294 * tau_b
    #Clear sky horizontal diffuse radiation
    I_cd = I_on * tau_d * cos(theta_z)
    #Total horizontal clear sky radiation
    I_th_cs = max(0, (I_cb + I_cd))
    #Simplified calculation (G.Masy)
    I_th_cs2 = max(0, (0.7 * I_on * cos(theta_z)))
    return I_th_cs,I_th_cs2

def Btest_cases(Btest=None, h=None):
    if (Btest == 195 or Btest == 200):
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_gl = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 210 or Btest == 220:
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 230:
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 1

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings
        #Ventilation
        V_dot_vent = 0

    elif Btest == 240:
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 250:
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        alpha_lopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 270 or Btest == 300):
        #set points
        T_i_set_h = 20
        T_i_set_c = 20

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 320:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 395:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_gl = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 400:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 410:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 0

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif Btest == 420:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = pd.Series([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings
        #Ventilation
        V_dot_vent = 0

    elif Btest == 430:
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 600 or Btest == 620 or Btest == 800 or Btest == 900 or Btest == 920):
        #set points
        T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 640 or Btest == 940):
        #set points
        if h > 23 or h <= 7:
            T_i_set_h = 10
        else:
            T_i_set_h = 20
        T_i_set_c = 27

        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = pd.Series([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = pd.Series([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = pd.Series([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = pd.Series([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = pd.Series([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = pd.Series([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 650 or Btest == 950):
        #set points
        if (h > 18 or h <= 7):
            V_dot_vent = 1703.16

    return T_i_set_h,T_i_set_c,SHGC_gl_0,epsilon_ir_hopw,epsilon_ir_lopw,epsilon_ir_gl,alpha_hopw,alpha_lopw,e_solshad,\
           mode_solshad,NL_ext_max,IAC_solshad, Q_dot_appl, ACH_inf, V_dot_vent

def DRYAIRPROP(t=None, P=None):
    #Dry air properties
    #Air gas constant
    R_a = 287# J/kg-K
    #Dry air compressibility factor
    Z_a = 0.97901976 + 1.8181227e-4 * (t + 273.15) - 5.2676556e-7 * (t + 273.15) ** 2 + \
          5.2528153e-10 * (t + 273.15) ** 3 - 6.6519834e-9 * P
    #dry air specific volume
    v_a = R_a * (t + 273.15) / P * Z_a
    # %dry air virial contribution to enthalpy
    # v_ha=1000*(3.1327784-3.1755794E-2*(t+273.15)+1.2006509E-4*(t+273.15)^2-1.5693167E-7*(t+273.15)^3-2.7847739E-6*P+6.951284E-14*P^2);
    # %Dry air enthalpy
    # h_a=1004.1923*t+v_ha;
    return v_a

def h_t_w(c_p_a=None, c_p_g=None, h_fg0=None, t=None, w=None):
    #Enthalpy as function of the drybulb temperature and the humidity ratio at
    #101325 Pa
    h = c_p_a * t + (h_fg0 + c_p_g * t) * w
    return h

def SOLARAZIMUTH(omega=None, delta=None, theta_z=None, phi=None):
    #SOLAR AZIMUTH ANGLE
    #Ref: Duffie,J.A.,Beckman, W.A. 1980. Solar engineering of thermal
    #processes. 2nd Edition. John Wiley & Sons.
    #p.16 & 17

    #INPUTS
    # -phi: Latitude in radians
    # -delta: Declination angle in radians
    # -omega: Hour angle in radians
    # -theta_z: Zenith angle in radians, i.e. angle of incidence of beam radiation on a horizontal surface

    #OUTPUT
    # -gamma_s: Solar azimuth angle in radians
    from math import sin,asin,tan,acos
    #pseudo azimuth angle
    singamma_s2 = (sin(omega) * sin(delta)) / (sin(theta_z))
    gamma_s2 = asin(singamma_s2)

    #Coefficients
    #Hour angle when the sun is due west or east
    cosomega_ew = tan(delta) / tan(phi + 1E-6)
    if abs(cosomega_ew) > 1:
        C1 = 1
    else:
        omega_ew = acos(cosomega_ew)
        if abs(omega) < omega_ew:
            C1 = 1
        else:
            C1 = -1

    if (phi * (phi - delta)) >= 0:
        C2 = 1
    else:
        C2 = -1

    if omega >= 0:
        C3 = 1
    else:
        C3 = -1
    gamma_s = C1 * C2 * gamma_s2 + C3 * 180 * ((1 - (C1 * C2))) / 2
    return gamma_s

def IR_horiz(I_ir_cs=None, I_ir_cc=None, J_cs=None, J_cc=None, I_glob_data=None, I_th_cs=None, tau_h=None,
             tau_h_1=None, tau_h_2=None):
#    try:
    if (I_glob_data[tau_h] > 0) and (I_th_cs[tau_h] > 0):
        J_rel = min(J_cs, max(J_cc, min(1, max(0, (I_glob_data[tau_h] / I_th_cs[tau_h])))))
    else:
        #Down value
        tau_down = tau_h
        while (I_th_cs[tau_down] == 0) and (I_glob_data[tau_down] == 0):
            if tau_down == tau_h_1:
                tau_down = tau_h_2
            tau_down = tau_down - 1
        #Up value
        tau_up = tau_h
        while (I_th_cs[tau_up] == 0) and (I_glob_data[tau_up] == 0):
            if tau_up == tau_h_2:
                tau_up = tau_h_1
            tau_up = tau_up + 1
        #Average value
        J_down = min(1, max(0, (I_glob_data[tau_down] / I_th_cs[tau_down])))
        J_up = min(1, max(0, (I_glob_data[tau_up] / I_th_cs[tau_up])))
        J_avg = (J_up + J_down) / 2
        J_rel = min(J_cs, max(J_cc, J_avg))
#    except:
#        print "I_ir_h error on timestep: "+ str(tau_up)
#        J_rel=0
#    Sky radiation (between 45 and 100 W/m)
    I_ir_h = -(I_ir_cc + (I_ir_cs - I_ir_cc) * (J_rel - J_cc) / (J_cs - J_cc))
    return I_ir_h

def INCIDENCEANG(phi=None, delta=None, omega=None, slope=None, surf_az=None):
    #ZENITH ANGLE
    #Ref: Duffie,J.A.,Beckman, W.A. 1980. Solar engineering of thermal
    #processes. 2nd Edition. John Wiley & Sons.
    #OUTPUTS
    # -gamma: Surface azimuth angle in radians
    # -beta: Surface slope in radians
    # -theta: Angle of incidence of beam radiation in radians
    #INPUTS
    # -phi: Latitude in radians
    # -delta: Declination angle in radians
    # -omega: Hour angle in radians
    # -surf_az: Surface azimuth -180<surf_az<180 (south: 0, east negative, west positive)
    # -slope: slope of the surface 0<slope<180
    from math import sin,cos,pi,acos
    #Angles in radians%
    gamma = surf_az * pi / 180
    beta = slope * pi / 180

    #Angle of incidence (between beam radiation on a surface and the normal of
    #the surface)%
    costheta1 = sin(delta) * sin(phi) * cos(beta)
    costheta2 = -sin(delta) * cos(phi) * sin(beta) * cos(gamma)
    costheta3 = cos(delta) * cos(phi) * cos(beta) * cos(omega)
    costheta4 = cos(delta) * sin(phi) * sin(beta) * cos(gamma) * cos(omega)
    costheta5 = cos(delta) * sin(beta) * sin(gamma) * sin(omega)
    theta = acos(costheta1 + costheta2 + costheta3 + costheta4 + costheta5)
    return theta

def DIRRAD(I_dir_n=None, theta=None, f_low_dir=None):
    #DIRECT RADIATION INCIDENT TO WALL
    #Ref: Duffie,J.A.,Beckman, W.A. 1980. Solar engineering of thermal
    #processes. 2nd Edition. John Wiley & Sons.

    #OUTPUTS
    # - I_dir_w: Direct solar radiation incident to wall

    #INPUTS
    # - theta: Direct solar radiation incidence angle (in radians)
    # - I_dir_n: Normal direct solar radiation
    # - f_low_dir: Direct solar radiation lowering factor (masks)

    #Direct solar radiation incident to the wall
    from math import cos
    I_dir_i = I_dir_n * max(0, cos(theta))
    #Beam radiation reaching the wall
    I_dir_w = f_low_dir * I_dir_i
    return I_dir_w

def HEMISRAD(albedo=None, I_glob_h=None, I_diff_h=None, I_dir_h=None, theta_z=None, slope=None, f_low_diff=None):
    #HEMISPHERICAL RADIATION INCIDENT TO WALL
    #Ref: Duffie,J.A.,Beckman, W.A. 1980. Solar engineering of thermal
    #processes. 2nd Edition. John Wiley & Sons.

    #OUTPUTS
    # - I_dir_h: Direct solar radiation on horizontal surface
    # - I_hemis_w: Hemispherical solar radiation incident to wall

    #INPUTS
    # - albedo: Ground albedo
    # - I_glob_h: Global solar radiation on horizontal surface
    # - I_diff_h
    # - theta_z: solar zenith angle (in radians)
    # - slope: Wall slope (in degrees)
    # - f_low_diff: diffuse radiation lowering factors (mask)
    from math import pi,cos
    #Angles in radians
    beta = slope * pi / 180

    #Horizontal diffuse solar radiation
    #I_dir_h=I_glob_h-I_diff_h;
    #Horizontal direct solar radiation
    #I_dir_n=I_dir_h/cos(theta_z);
    #Diffuse solar radiation incident to the wall
    I_diff_i = 1 / 2 * I_diff_h * (1 + cos(beta))
    #Reflected solar radiation incident to the wall
    I_ref_i = (I_dir_h + I_diff_h) * 1 / 2 * albedo * (1 - cos(beta))
    #Diffuse radiation reaching the wall
    I_diff_w = f_low_diff * I_diff_i
    #Reflected radiation reaching the wall
    I_ref_w = I_ref_i
    #Hemispherical radiation reaching the wall
    I_hemis_w = I_diff_w + I_ref_w
    return I_ref_w, I_diff_w, I_hemis_w

def INCIDENTRAD(Lat=None, Long=None, Long_st=None, albedo=None, n=None, h=None, I_glob_h=None, I_diff_h=None,
                surf_az=None, slope=None, f_low_dir=None, f_low_diff=None):
    #INCIDENT RADIATION FOR GIVEN ORIENTATION AND SLOPELat,Long,Long_st,albedo,n,h,I_glob_h,I_diff_h,surf_az(ori),slope(ori),f_low_diff(ori),f_low_dir(ori)
    #OUTPUTS
    # - theta: Incidence angle of direct radiation for given orientation/slope (rad)
    # - I_tot_w: Total solar radiation incident to wall (W/m)
    # - I_ref_w: Reflected solar radiation incident to wall (W/m)
    # - I_diff_w: Diffuse solar radiation incident to wall (W/m)
    # - I_dir_w: Direct solar radiation incident to wall (W/m)

    #INPUTS
    # - Lat: Latitude of the location (north positive) -90<Lat<90
    # - Long: Longitude of the location (west positive) 0<Long<180
    # - Long_st: Longitude of the standard meridian of the time zone
    # - albedo: Ground albdo
    # - n: day 1<n<365
    # - h: hour 1<h<8760
    # - I_glob_h: Global horizontal solar radiation
    # - I_dir_n: Normal direct solar radiation
    # - surf_az: Surface azimuth angle for given orientation
    # - slope: Surface slope
    from math import cos
    if I_glob_h > 0:
        #Main angles for location
        phi, delta, omega, theta_z, h_sol = ZENITHANG(Lat, Long, Long_st, n, h)
        #Solar azimuth angle
        gamma_s = SOLARAZIMUTH(omega, delta, theta_z, phi)
        #Direct solar radiation on horizontal surface
        I_dir_h = max(0, I_glob_h - I_diff_h)
        #Direct normal solar radiation
        if cos(theta_z) > 0.0523:
            I_dir_n = max(0, I_dir_h / max(0, min(1, cos(theta_z))))
        else:
            I_dir_n = 0
        #DIRECT AND HEMISPHERICAL RADIATIONS
        #Incidence angle of direct radiation for given orientation/slope(radians)
        theta = INCIDENCEANG(phi, delta, omega, slope, surf_az)
        #Incident direct solar radiation for given orientation/slope(W/m)
        I_dir_w = DIRRAD(I_dir_n, theta, f_low_dir)
        #Incident hemispherical solar radiation for given orientation/slope(W/m)
        I_ref_w, I_diff_w, I_hemis_w = HEMISRAD(albedo, I_glob_h, I_diff_h, I_dir_h, theta_z, slope, f_low_diff)
        #Total incident radiation for main orientations
        I_tot_w = I_ref_w + I_diff_w + I_dir_w

    else:
        theta_z = 0
        theta = 0
        gamma_s = 0
        I_tot_w = 0
        I_ref_w = 0
        I_diff_w = 0
        I_dir_w = 0

    return theta, theta_z, gamma_s, I_tot_w, I_ref_w, I_diff_w, I_dir_w, I_dir_n

def SOLARSHAD(mode_solshad=None, NL_ext_max=None, I_tot=None):
    #FACTOR OF USE OF SOLAR SHADINGS
    #Reference:
    # - Th-CE 2005 Version 7.3. 15/03/2006
    # - ALESSANDRINI J.M., FLEURY, E, FILFLI S., MARCHIO D. 2006. Impact de la gestion de lclairage et des protections solaires sur la
    # consommation dnergie de btiments de bureaux climatiss.Climamed,Lyon,France.

    #INPUTS
    # - I_tot: total solar radiation incident to the wall (W/m)

    #OUTPUT
    # - f_use: factor of use of the external solar shadings

    NL_ext = I_tot * 100#100 lm/W

    if mode_solshad == 1:
        #Manual solar shadings
        if NL_ext < 25000:
            f_shad = 0.07 + (0.25 - 0.07) / (25000 - 0) * (NL_ext - 0)
        elif (NL_ext >= 25000) and (NL_ext < 100000):
            f_shad = 0.25 + (0.45 - 0.25) / (100000 - 25000) * (NL_ext - 25000)
        else:
            f_shad = 0.45
    else:
        #Automatic solar shadings
        if NL_ext > NL_ext_max:
            f_shad = 1
        else:
            f_shad = 1 / NL_ext_max * NL_ext
    return f_shad, NL_ext

def VARSHGC(theta=None, p=None):
    from math import pi,tan
    zero = 1E-6
    if (theta > pi / 2) and (theta < 3 * pi / 2):
        k_SHGC = 0
    else:
        if (p < zero):
            k_SHGC = 1
        else:
            k_SHGC = min(1, (1 - (tan(theta / 2)) ** p))
    return k_SHGC

def SOLARGAINS(Lat=None, Long=None, Long_st=None, albedo=None, n_walls=None, ori=None, SHGC_gl_0=None, p_SHGC=None,
               f_sg11=None, f_sg12=None, f_sg21=None, f_sg22=None, f_hemis=None, A_v_hopw=None, A_v_lopw=None,
               A_t_gl=None, surf_az=None, slope=None, f_low_dir=None, f_low_diff=None, zero=None, alpha_hopw=None,
               alpha_lopw=None, e_solshad=None, mode_solshad=None, IAC_solshad=None, NL_ext_max=None, n=None, h=None,
               I_glob_h=None, I_diff_h=None, I_ir_h=None, epsilon_ir_lopw=None, epsilon_ir_gl=None,
               epsilon_ir_hopw=None, U_lopw=None, U_gl=None, U_hopw=None, h_e_l=None, h_e_h=None):
    #OUTPUTS

    # - I_sol_rf: Solar radiation absorbed by horizontal roof(W/m)
    # - I_sol_lopw: Solar radiation absorbed by vertical light opaque walls(W/m)
    # - I_sol_hopw: Solar radiation absorbed by vertical heavy opaque walls(W/m)
    # - I_sol_gl: Solar radiation transmitted through vertical windows(W/m)
    # - Q_dot_sl: total solar flux and infrared emission of light opaque walls(W)
    # - Q_dot_sh: total solar flux and infrared emission of heavy opaque walls
    # - Q_dot_sd_tot: total direct solar flux through glazing
    # - Q_dot_svl_tot: Solar gain due to the increase in temperature of the ventilated air cavity


    #INPUTS
    # - Lat: Latitude of the location (north positive) -90<Lat<90
    # - Long: Longitude of the location (west positive) 0<Long<180
    # - Long_st: Longitude of the standard meridian of the time zone
    # - albedo: Ground albdo
    # - SHGC_gl_0: Normal solar heat gain coefficient of window
    # - p_SHGC: Angular dependency factor of SHGC for direct radiation
    # - A_v_hopw,A_v_lopw: Vertical Opaque Walls surface area for each orientation (m)
    # - A_t_gl: Vertical Windows surface area for each orientation (m)
    # - n: day 1<n<365
    # - h: hour 1<h<8760
    # - I_glob_h: Global horizontal solar radiation
    # - I_diff_h: Diffuse horizontal solar radiation

    #Initialization:
    # theta=zeros(1, n);
    # theta_z=zeros(1, n_ori);
    # gamma_s=zeros(1, n_ori);
    # I_tot_w=zeros(1, n_ori);
    # I_diff_w=zeros(1, n_ori);
    # I_ref_w=zeros(1, n_ori);
    # I_dir_w=zeros(1, n_ori);
    # I_dir_n=zeros(1, n_ori);
    # I_hemis_w=zeros(1, n_ori);
    # I_sol_lopw=zeros(1, n_ori);
    # I_sol_hopw=zeros(1, n_ori);
    # Q_dot_sol_l=zeros(1, n_ori);
    # Q_dot_sol_h=zeros(1, n_ori);
    from numpy import zeros
    from math import pi,cos
    import pandas as pd

    #OLD ATTEMPT TO CONTAIN SURFACE CALCULATIONS IN LISTS - HAD MESSY INITIALIZATION PROBLEMS
    #    Q_dot_sol_h, Q_dot_sol_l, Q_dot_sol_gl_dir, Q_dot_sol_gl_hemis, Q_dot_sd_dir, Q_dot_sd_hemis, Q_dot_svl_dir,\
    #    Q_dot_svl_hemis, I_tot_w, I_hemis_w = [],[],[],[],[],[],[],[],[],[]
    #
    #    theta, theta_z, gamma_s, I_tot_w, I_ref_w, I_diff_w, I_dir_w, I_dir_n = [],[],[],[],[],[],[],[]
    #
    #    f_shad, NL_ext, I_sol_hopw, I_sol_lopw, SHGC1, SHGC2, SHGC_dir1, SHGC_hemis1, SHGC_dir2, SHGC_hemis2 \
    #    = [],[],[],[],[],[],[],[],[],[]
    #
    #    Sp11_dir, Sp21_dir, Sp31_dir, Sp11_hemis, Sp21_hemis, Sp31_hemis, Sp12_dir, Sp22_dir, SHGC_dir2, Sp32_dir, \
    #    SHGC_dir2, Sp12_dir, Sp22_dir, Sp12_hemis, SHGC_hemis2, Sp22_hemis, Sp32_hemis = [],[],[],[],[],[],[],[],[],[],[],\
    #    [],[],[],[],[],[]
    #
    #    Q_dot_IR_l, Q_dot_IR_h = {},{}

    SurfaceDataFrameLabels = ['Q_dot_sol_h',
                              'I_diff_w',
                              'Q_dot_sol_gl_dir',
                              'SHGC_dir1',
                              'SHGC_dir2',
                              'Sp32_hemis',
                              'Sp31_dir',
                              'Sp22_hemis',
                              'Q_dot_sol_l',
                              'Q_dot_sd_hemis',
                              'Q_dot_svl_hemis',
                              'I_sol_hopw',
                              'I_dir_w',
                              'NL_ext',
                              'Sp11_dir',
                              'Sp21_dir',
                              'Sp22_dir',
                              'Sp12_dir',
                              'theta',
                              'I_dir_n',
                              'gamma_s',
                              'I_tot_w',
                              'I_hemis_w',
                              'SHGC_hemis2',
                              'Q_dot_svl_dir',
                              'SHGC1',
                              'SHGC2',
                              'Q_dot_sd_dir',
                              'Q_dot_IR_l',
                              'I_sol_lopw',
                              'Sp11_hemis',
                              'f_shad',
                              'Sp21_hemis',
                              'Sp31_hemis',
                              'Q_dot_sol_gl_hemis',
                              'SHGC_hemis1',
                              'Sp12_hemis',
                              'Sp32_dir',
                              'Q_dot_IR_h',
                              'I_ref_w',
                              'theta_z']

    SurfaceList = [(x) for x in range(0,n_walls)]
    Surface = pd.DataFrame(columns=SurfaceDataFrameLabels,index=SurfaceList)

    for i in range(0,n_walls):
        if I_glob_h > 0:
            if (ori[i] == 'F') or (ori[i] == 'C'):
                Surface.Q_dot_sol_h[i] = 0
                Surface.Q_dot_sol_l[i] = 0
                Surface.Q_dot_sol_gl_dir[i] = 0
                Surface.Q_dot_sol_gl_hemis[i] = 0
                Surface.Q_dot_sd_dir[i] = 0
                Surface.Q_dot_sd_hemis[i] = 0
                Surface.Q_dot_svl_dir[i] = 0
                Surface.Q_dot_svl_hemis[i] = 0
                Surface.I_tot_w[i] = 0
            else:
                Surface.theta[i], Surface.theta_z[i], Surface.gamma_s[i], Surface.I_tot_w[i], Surface.I_ref_w[i],\
                Surface.I_diff_w[i], Surface.I_dir_w[i], Surface.I_dir_n[i] = \
                INCIDENTRAD(Lat, Long, Long_st, albedo, n, h, I_glob_h, I_diff_h, surf_az[i], slope[i], f_low_dir[i],f_low_diff[i])

                Surface.I_hemis_w[i] = Surface.I_ref_w[i] + Surface.I_diff_w[i]

                #USE OF SHADINGS - Fraction of the area covered by shading system
                #(vertical openings only)
                if e_solshad[i] == 0:
                    Surface.f_shad[i] = 0
                    Surface.NL_ext[i] = 0
                else:
                    Surface.f_shad[i], Surface.NL_ext[i] = SOLARSHAD(mode_solshad[i], NL_ext_max[i], Surface.I_tot_w[i])

                # VERTICAL WALLS AND ROOF
                #Heavy opaque vertical walls incident radiation
                Surface.I_sol_hopw[i] = Surface.I_tot_w[i]            #(W)
                Surface.Q_dot_sol_h[i] = alpha_hopw[i] * A_v_hopw[i] * Surface.I_sol_hopw[i]
                #Light opaque walls incident radiation
                Surface.I_sol_lopw[i] = Surface.I_tot_w[i]            #(W)
                Surface.Q_dot_sol_l[i] = alpha_lopw[i] * A_v_lopw[i] * Surface.I_sol_lopw[i]

                #VERTICAL AND HORIZONTAL WINDOWS
                #SHGC and Light Transmittance without exterior solar shading
                Surface.SHGC1[i] = SHGC_gl_0[i]
                #SHGC and Light Transmittance with exterior solar shading
                Surface.SHGC2[i] = IAC_solshad[i] * SHGC_gl_0[i]

                #Glazed walls incident radiation
                #SHGC for direct solar radiation without solar shadings
                Surface.SHGC_dir1[i] = Surface.SHGC1[i] * VARSHGC(Surface.theta[i], p_SHGC)
                #SHGC for hemispherical solar radiation without solar shadings
                Surface.SHGC_hemis1[i] = Surface.SHGC1[i] * f_hemis
                #SHGC for direct solar radiation with solar shadings
                Surface.SHGC_dir2[i] = Surface.SHGC2[i] * VARSHGC(Surface.theta[i], p_SHGC)
                #SHGC for hemispherical solar radiation with solar shadings
                Surface.SHGC_hemis2[i] = Surface.SHGC2[i] * f_hemis

                #Components of the solar factor of a window
                # without solar protection
                Surface.Sp11_dir[i] = f_sg11 * Surface.SHGC_dir1[i]       #Assumption when no information on the glazing
                Surface.Sp21_dir[i] = f_sg21 * Surface.SHGC_dir1[i]
                Surface.Sp31_dir[i] = Surface.SHGC_dir1[i] - Surface.Sp11_dir[i] - Surface.Sp21_dir[i]
                Surface.Sp11_hemis[i] = f_sg11 * Surface.SHGC_hemis1[i]   #Assumption when no information on the glazing
                Surface.Sp21_hemis[i] = f_sg21 * Surface.SHGC_hemis1[i]
                Surface.Sp31_hemis[i] = Surface.SHGC_hemis1[i] - Surface.Sp11_hemis[i] - Surface.Sp21_hemis[i]
                # with solar ventilated protection
                Surface.Sp12_dir[i] = f_sg12 * Surface.SHGC_dir2[i]
                Surface.Sp22_dir[i] = f_sg22 * Surface.SHGC_dir2[i]
                Surface.Sp32_dir[i] = Surface.SHGC_dir2[i] - Surface.Sp12_dir[i] - Surface.Sp22_dir[i]
                Surface.Sp12_hemis[i] = f_sg12 * Surface.SHGC_hemis2[i]
                Surface.Sp22_hemis[i] = f_sg22 * Surface.SHGC_hemis2[i]
                Surface.Sp32_hemis[i] = Surface.SHGC_hemis2[i] - Surface.Sp12_hemis[i] - Surface.Sp22_hemis[i]

                Surface.Q_dot_sol_gl_dir[i] = A_t_gl[i] * Surface.I_dir_w[i] * ((1 - Surface.f_shad[i]) *
                                                                                Surface.Sp21_dir[i] + Surface.f_shad[i]
                                                                                * Surface.Sp22_dir[i])
                Surface.Q_dot_sol_gl_hemis[i] = A_t_gl[i] * Surface.I_hemis_w[i] * ((1 - Surface.f_shad[i]) *
                                                                                    Surface.Sp21_hemis[i] + Surface.f_shad[i]
                                                                                    * Surface.Sp22_hemis[i])
                #Direct solar gain through windows
                Surface.Q_dot_sd_dir[i] = Surface.I_dir_w[i] * A_t_gl[i] * ((1 - Surface.f_shad[i]) * Surface.Sp11_dir[i]
                                                                            + Surface.f_shad[i] * Surface.Sp12_dir[i])
                Surface.Q_dot_sd_hemis[i] = Surface.I_hemis_w[i] * A_t_gl[i] * ((1 - Surface.f_shad[i]) * Surface.Sp11_hemis[i]
                                                                                + Surface.f_shad[i] * Surface.Sp12_hemis[i])
                #Solar gain due to the increase in temperature of the
                #ventilated air cavity
                Surface.Q_dot_svl_dir[i] = A_t_gl[i] * (Surface.Sp31_dir[i] * (1 - Surface.f_shad[i]) * Surface.I_dir_w[i]
                                                        + Surface.Sp32_dir[i] * Surface.f_shad[i] * Surface.I_dir_w[i])
                Surface.Q_dot_svl_hemis[i] = A_t_gl[i] * (Surface.Sp31_hemis[i] * (1 - Surface.f_shad[i]) * Surface.I_hemis_w[i]
                                                          + Surface.Sp32_hemis[i] * Surface.f_shad[i] * Surface.I_hemis_w[i])

                Q_dot_sol_h_tot = Surface.Q_dot_sol_h.sum()
                Q_dot_sol_l_tot = Surface.Q_dot_sol_l.sum()
                Q_dot_sol_gl_tot = Surface.Q_dot_sol_gl_dir.sum() + Surface.Q_dot_sol_gl_hemis.sum()
                Q_dot_sd_tot = Surface.Q_dot_sd_dir.sum() + Surface.Q_dot_sd_hemis.sum()
                Q_dot_svl_tot = Surface.Q_dot_svl_dir.sum() + Surface.Q_dot_svl_hemis.sum()

        else:
#            Surface.I_tot_w = zeros((1, n_walls))
            Surface.I_tot_w = pd.Series([0 for x in range(0,n_walls)])
            Q_dot_sol_h_tot = 0
            Q_dot_sol_l_tot = 0
            Q_dot_sol_gl_tot = 0
            Q_dot_sd_tot = 0
            Q_dot_svl_tot = 0

        #INFRARED EMISSION
        if (ori[i]=='F') or (ori[i]=='C'):
            #light walls
            #Q_dot_IR_l(ori)=(A_v_lopw(ori)*1/2*(1+cos(beta))*epsilon_ir_lopw*I_ir_h*U_lopw+A_t_gl(ori)*1/2*(1+cos(beta))*epsilon_ir_gl*I_ir_h*U_gl)/h_e_l ;
            # Heavy opaque walls
            #Q_dot_IR_h(ori)=A_v_hopw(ori)*epsilon_ir_hopw*I_ir_h*1/2*(1+cos(beta))*U_hopw/h_e_h;
            #light walls
            Surface.Q_dot_IR_l[i] = 0
            # Heavy opaque walls
            Surface.Q_dot_IR_h[i] = 0
        else:
            beta_w = slope[i] * pi / 180
            #light walls
            Surface.Q_dot_IR_l[i] = (A_v_lopw[i] * 1 / 2 * (1 + cos(beta_w)) * epsilon_ir_lopw[i] * I_ir_h +
                                    A_t_gl[i] * 1/2 * (1 + cos(beta_w)) * epsilon_ir_gl[i] * I_ir_h)
            # Heavy opaque walls
            Surface.Q_dot_IR_h[i] = A_v_hopw[i] * epsilon_ir_hopw[i] * I_ir_h * 1 / 2 * (1 + cos(beta_w))
        
    Q_dot_IR_l_tot = Surface.Q_dot_IR_l.sum()
    Q_dot_IR_h_tot = Surface.Q_dot_IR_h.sum()

    #SOLAR GAINS
    Q_dot_sl_tot = Q_dot_sol_l_tot + Q_dot_sol_gl_tot - Q_dot_IR_l_tot
    Q_dot_sh_tot = Q_dot_sol_h_tot - Q_dot_IR_h_tot

    #Total transmitted solar radiation
    I_tr_tot = (Q_dot_sol_gl_tot + Q_dot_sd_tot + Q_dot_svl_tot) / (sum(A_t_gl) + zero)

    return Q_dot_sl_tot, Q_dot_sh_tot, Q_dot_sd_tot, Q_dot_svl_tot, Surface.I_tot_w, I_tr_tot, Q_dot_IR_l_tot, Q_dot_IR_h_tot, \
           Q_dot_sol_l_tot, Q_dot_sol_h_tot, Q_dot_sol_gl_tot

def TEMP_OUT(error=None, zero=None, Q_dot_sl=None, Q_dot_sh=None, T_e=None, H_tr_em=None, H_tr_es=None, H_tr_ms=None,
             h_e_l=None, h_e_h=None, A_lopw_t=None, A_gl_t=None, A_fr_t=None, A_hopw_t=None):
    if abs(Q_dot_sl) < error:
        T_es = T_e
    else:
        T_es = T_e + Q_dot_sl / (h_e_l * (A_lopw_t + A_gl_t + A_fr_t) + zero)
        #T_es=T_e+Q_dot_sl/(H_tr_es+zero);

    if abs(Q_dot_sh) < error:
        T_em = T_e
    else:
        #T_em=T_e+Q_dot_sh/(H_tr_em+H_tr_ms+zero);
        T_em = T_e + Q_dot_sh / (h_e_h * A_hopw_t + zero)
    return T_es, T_em
