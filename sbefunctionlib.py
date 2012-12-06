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
    
    return (n_walls, f_low_diff, f_low_dir, ori, surf_az, slope, A_t, A_fl, A_lopw_t, A_hopw_t, A_gl_t, A_fr_t, A_lopw,
            A_hopw, A_gl, h_cl, C_m, A_m, U_hopw, U_lopw, U_fr, U_gl)

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_gl = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings
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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        alpha_lopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        epsilon_ir_gl = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        alpha_lopw = mcat([0.1, 0.1, 0.1, 0.1, 0.1, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings
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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

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
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 640 or Btest == 940):
        #set points
        if h > 23 or h <= 7:
            T_i_set_h = 10
        else:
            T_i_set_h = 20
        end
        T_i_set_c = 27


        #internal gain
        Q_dot_appl = 200

        #infiltrations
        ACH_inf = 0.5

        #SOLAR PROPERTIES
        #SHGC
        SHGC_gl_0 = mcat([0.789, 0.789, 0.789, 0.789, 0.789, 0, 0])
        #IR emittance
        epsilon_ir_hopw = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        epsilon_ir_lopw = mcat([0, 0, 0, 0, 0, 0, 0])
        epsilon_ir_gl = mcat([0.9, 0.9, 0.9, 0.9, 0.9, 0, 0])
        #Solar absorbance
        alpha_hopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        alpha_lopw = mcat([0.6, 0.6, 0.6, 0.6, 0.6, 0, 0])
        #Solar Shadings
        e_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #0=no solar shading; 1=interior solar shadings; 2=exterior solar shadings
        mode_solshad = mcat([1, 1, 1, 1, 1, 0, 0])    #1=manual solar shadings; 2=automatic solar shadings
        NL_ext_max = mcat([0, 0, 0, 0, 0, 0, 0])    #Exterior natural lighting intensity for control of shadings
        IAC_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Indoor solar Attenuation Coefficient (fraction of SHGC with solar shadings)
        f_c_solshad = mcat([0, 0, 0, 0, 0, 0, 0])    #Convective fraction of solar gains with solar shadings

        #Ventilation
        V_dot_vent = 0

    elif (Btest == 650 or Btest == 950):
        #set points
        if (h > 18 or h <= 7):
            V_dot_vent = 1703.16
