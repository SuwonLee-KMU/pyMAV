import numpy as np
from math import cos, sin

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

class airframe(object):
    def __init__(self, name):
        self.name = name # 'zagi', 'aerosonde'
        if name == 'zagi':
            self.zagi()
        elif name == 'aerosonde':
            self.aerosonde()
        else:
            raise ValueError("Invalid name of the airframe.")
    
    def zagi(self):
        self.MASS = 1.56 # kg
        self.Jx = 0.1147 # kg m^2
        self.Jy = 0.0576 # kg m^2
        self.Jz = 0.1712 # kg m^2
        self.Jxz = 0.0015 # kg m^2
        self.S = 0.2589 # m^2
        self.b = 1.4224 # m
        self.c = 0.3302 # m
        self.S_prop = 0.0314 # m^2
        self.C_prop = 1.0
        self.rho = 1.2682 # kg/m^3 <- why here...?
        self.k_motor = 20
        self.k_T_p = 0
        self.k_omega = 0
        self.e = 0.9
        self.CL0 = 0.09167
        self.CD0 = 0.01631
        self.Cm0 = -0.02338
        self.CLa = 3.5016
        self.CDa = 0.2108
        self.Cma = -0.5675
        self.CLq = 2.8932
        self.CDq = 0
        self.Cmq = -1.3990
        self.CLde = 0.2724
        self.CDde = 0.3045
        self.Cmde = -0.3254
        self.CDp = 0.0254
        self.CY0 = 0
        self.Cl0 = 0
        self.Cn0 = 0
        self.CYb = -0.07359
        self.Clb = -0.02854
        self.Cnb = -0.00040
        self.CYp = 0
        self.Clp = -0.3209
        self.Cnp = -0.01297
        self.CYr = 0
        self.Clr = 0.03066
        self.Cnr = -0.00434
        self.CYda = 0
        self.Clda = 0.1682
        self.Cnda = -0.00328

        # What are you...?
        self.M = 50
        self.alpha0 = 0.4712
        self.epsilon = 0.1592
    
    def aerosonde(self):
        self.MASS = 13.5 # kg
        self.Jx = 0.8244 # kg m^2
        self.Jy = 1.135 # kg m^2
        self.Jz = 1.759 # kg m^2
        self.Jxz = 0.1204 # kg m^2
        self.S = 0.55 # m^2
        self.b = 2.8956 # m
        self.c = 0.18994 # m
        self.S_prop = 0.2027 # m^2
        self.rho = 1.2682 # kg/m^3
        self.k_motor = 80
        self.k_T_p = 0
        self.k_omega = 0
        self.e = 0.9
        self.CL0 = 0.28
        self.CD0 = 0.03
        self.Cm0 = -0.02338
        self.CLa = 3.45
        self.CDa = 0.30
        self.Cma = -0.38
        self.CLq = 0
        self.CDq = 0
        self.Cmq = -3.6
        self.CLde = -0.36
        self.CDde = 0
        self.Cmde = -0.5
        self.C_prop = 1.0
        self.CDp = 0.0437
        self.Cndr = -0.032
        self.CY0 = 0
        self.Cl0 = 0
        self.Cn0 = 0
        self.CYb = -0.98
        self.Clb = -0.12
        self.Cnb = 0.25
        self.CYp = 0
        self.Clp = -0.26
        self.Cnp = 0.022
        self.CYr = 0
        self.Clr = 0.14
        self.Cnr = -0.35
        self.CYda = 0
        self.Clda = 0.08
        self.CYdr = -0.17
        self.Cldr = 0.105
        self.Cnda = 0
        
        self.M = 50
        self.alpha0 = 0.4712
        self.epsilon = 0.1592

    def get_aero_coef(self):
        ac = dotdict()
        attr = self.__dict__
        for i, v in enumerate(attr):
            if v.startswith('C'):
                ac[v] = getattr(self, v)
        return ac

    def get_airframe_param(self):
        ap = dotdict()
        attr = self.__dict__
        for i, v in enumerate(attr):
            if not v.startswith('C') :
                ap[v] = getattr(self, v)
        return ap

    def get_inertia(self):
        MOMENT_OF_INERTIA = [self.Jx, self.Jy, self.Jz]
        PRODUCT_OF_INERTIA = [0, self.Jxz, 0]
        return MOMENT_OF_INERTIA, PRODUCT_OF_INERTIA
    
    def get_lift_coefficients(self):
        return [self.CL0, self.CLa, self.CLq, self.CLde]
    
    def get_drag_coefficients(self):
        return [self.CD0, self.CDa, self.CDq, self.CDde]

    def get_pitching_moment_coefficients(self):
        return [self.Cm0, self.Cma, self.Cmq, self.Cmde]

    def get_lateral_force_coefficients(self):
        return [self.CY0, self.CYb, self. CYp, self.CYr, self.CYda, self.CYdr]

    def get_roll_moment_coefficients(self):
        return [self.Cl0, self.Clb, self.Clp, self.Clr, self.Clda, self.Cldr]

    def get_yaw_moment_coefficients(self):
        return [self.Cn0, self.Cnb, self.Cnp, self.Cnr, self.Cnda, self.Cndr]

    def get_inertia_matrix(self):
        MOMENT_OF_INERTIA, PRODUCTS_OF_INERTIA = self.get_inertia()
        Jx, Jy, Jz = np.reshape(MOMENT_OF_INERTIA, (3))
        Jxy, Jxz, Jyz = np.reshape(PRODUCTS_OF_INERTIA, (3))
        return [[Jx, -Jxy, -Jxz],
                [-Jxy, Jy, -Jyz],
                [-Jxz, -Jyz, Jz]]

    def get_Gamma_components(self):
        MOMENT_OF_INERTIA, PRODUCTS_OF_INERTIA = self.get_inertia()
        Jx, Jy, Jz = np.reshape(MOMENT_OF_INERTIA, (3))
        Jxy, Jxz, Jyz = np.reshape(PRODUCTS_OF_INERTIA, (3))
        Gamma = Jx*Jz - Jxz*Jxz
        Gamma1 = Jxz*(Jx - Jy + Jz) / Gamma
        Gamma2 = (Jz*(Jz-Jy) + Jxz*Jxz) / Gamma
        Gamma3 = Jz / Gamma
        Gamma4 = Jxz / Gamma
        Gamma5 = (Jz - Jx) / Jy
        Gamma6 = Jxz / Jy
        Gamma7 = ((Jx-Jy)*Jx + Jxz*Jxz) / Gamma
        Gamma8 = Jx / Gamma
        return Gamma, (Gamma1, Gamma2, Gamma3, Gamma4, Gamma5, Gamma6, Gamma7, Gamma8)

    def get_linear_coef(self):
        lon_linear_coef = [
            self.get_lift_coefficients(),
            self.get_drag_coefficients(),
            self.get_pitching_moment_coefficients()
        ]
        lat_linear_coef = [
            self.get_lateral_force_coefficients(),
            self.get_roll_moment_coefficients(),
            self.get_yaw_moment_coefficients()
        ]
        return lon_linear_coef, lat_linear_coef

    def get_pr_coefficients(self):
        _, Gamma_i = self.get_Gamma_components()
        Cl0, Clb, Clp, Clr, Clda, Cldr = self.get_roll_moment_coefficients()
        Cn0, Cnb, Cnp, Cnr, Cnda, Cndr = self.get_yaw_moment_coefficients()
        G1, G2, G3, G4, G5, G6, G7, G8 = Gamma_i
        Cp0 = G3*Cl0 + G4*Cn0
        Cpb = G3*Clb + G4*Cnb
        Cpp = G3*Clp + G4*Cnp
        Cpr = G3*Clr + G4*Cnr
        Cpda = G3*Clda + G4*Cnda
        Cpdr = G3*Cldr + G4*Cldr
        Cr0 = G4*Cl0 + G8*Cn0
        Crb = G4*Clb + G8*Cnb
        Crp = G4*Clp + G8*Cnp
        Crr = G4*Clr + G8*Cnr
        Crda = G4*Clda + G8*Cnda
        Crdr = G4*Cldr + G8*Cndr
        Cp = [Cp0, Cpb, Cpp, Cpr, Cpda, Cpdr]
        Cr = [Cr0, Crb, Crp, Crr, Crda, Crdr]
        return Cp, Cr
    
    def compute_XZ_coefficients(self, AOA):
        CL = self.get_lift_coefficients()
        CD = self.get_drag_coefficients()
        CX0, CXa, CXq, CXde = [-i*cos(AOA) + j*sin(AOA) for i, j in zip(CD, CL)]
        CZ0, CZa, CZq, CZde = [-i*sin(AOA) - j*cos(AOA) for i, j in zip(CD, CL)]
        CX = [CX0, CXa, CXq, CXde]
        CZ = [CZ0, CZa, CZq, CZde]
        return CX, CZ

class dynamics(object):
    pass