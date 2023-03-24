
def solve_Vsa(Vsa, Vsa0, Csa, Csa_l, Rs_u, Rs_l, Csa_u, Pext_l, Psa_u, Psv_u, 
              Psv_l, rho, g, H, Q):
    """ODE to solve for systemic arterial volume. Arguments are named below
    For more detailed information about units see the parameters.py file

    Args:
        Vsa (float): systemic arterial volume
        Vsa0 (float): reserve arterial volume
        Csa (float): systemic arterial compliance
        Csa_l (float): systemic upper arterial compliance
        Rs_u (float): systemic upper venous resistance
        Rs_l (float): systemic lower venous resistance
        Csa_u (float): systemic upper arterial compliance
        Pext_l (float): lower pressure external to vascular system
        Psa_u (float): pressure in upper systemic arteries
        Psv_u (float): pressure in upper systemic veins
        Psv_l (float): pressure in lower systemic veins
        rho (float): density of blood
        g (float): value of gravity
        H (float): total height
        Q (float): cardiac output

    Returns:
        float: dVsa/dt[t]
    """
    dydt = Q - Vsa/Csa*(1/Rs_u + 1/Rs_l) - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_l*rho*g*H)/Csa) - Psv_u) / Rs_u - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_u*rho*g*H)/Csa
            ) - Psv_l) / Rs_l # dVsa/dt
    
    return dydt

def solve_Vsv(Vsa, Vsa0, Csa, Csa_l, Rs_u, Rs_l, Csa_u, Pext_l, Psa_u, Psv_u, 
              Psv_l, rho, g, H, Q):
    """ODE to solve for systemic venous volume. Arguments are named below
    For more detailed information about units see the parameters.py file
    dVsv/dt = -dVsa/dt

    Args:
        Vsa (float): systemic arterial volume
        Vsa0 (float): reserve arterial volume
        Csa (float): systemic arterial compliance
        Csa_l (float): systemic upper arterial compliance
        Rs_u (float): systemic upper venous resistance
        Rs_l (float): systemic lower venous resistance
        Csa_u (float): systemic upper arterial compliance
        Pext_l (float): lower pressure external to vascular system
        Psa_u (float): pressure in upper systemic arteries
        Psv_u (float): pressure in upper systemic veins
        Psv_l (float): pressure in lower systemic veins
        rho (float): density of blood
        g (float): value of gravity
        H (float): total height
        Q (float): cardiac output

    Returns:
        float: dVsv/dt[t]
    """
    dydt = Q - Vsa/Csa*(1/Rs_u + 1/Rs_l) - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_l*rho*g*H)/Csa) - Psv_u) / Rs_u - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_u*rho*g*H)/Csa
            ) - Psv_l) / Rs_l # dVsv/dt
    return - dydt


def solve_ivp_Vsa(t, Vsa, Vsa0, Csa, Csa_l, Rs_u, Rs_l, Csa_u, Pext_l, Psa_u, Psv_u, 
              Psv_l, rho, g, H, Q):
    """ODE to solve for systemic arterial volume. Arguments are named below
    For more detailed information about units see the parameters.py file

    Args:
        Vsa (np.array): systemic arterial volume
        Vsa0 (float): reserve arterial volume
        Csa (float): systemic arterial compliance
        Csa_l (float): systemic upper arterial compliance
        Rs_u (float): systemic upper venous resistance
        Rs_l (float): systemic lower venous resistance
        Csa_u (float): systemic upper arterial compliance
        Pext_l (float): lower pressure external to vascular system
        Psa_u (float): pressure in upper systemic arteries
        Psv_u (float): pressure in upper systemic veins
        Psv_l (float): pressure in lower systemic veins
        rho (float): density of blood
        g (float): value of gravity
        H (float): total height
        Q (float): cardiac output

    Returns:
        float: dVsa/dt[t]
    """
    dydt = Q - Vsa[t]/Csa*(1/Rs_u + 1/Rs_l) - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_l*rho*g*H)/Csa) - Psv_u) / Rs_u - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_u*rho*g*H)/Csa
            ) - Psv_l) / Rs_l # dVsa/dt
    
    return dydt
