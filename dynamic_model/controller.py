import numpy as np

def get_linear_heart_rate(dPsa, F_max, F_min, P_max, P_min):
    """return a linear function of heart rate F in terms of pressure

    Args:
        F_max (_type_): maximum viable pressure for individual
        F_min (_type_): minimum viable pressure 
        P_max (_type_): maximum viable upper (coronary) 
            arterial pressure for individual
        P_min (_type_): minimum viable upper (coronary) 
            arterial pressure for individual
    """
    m = (F_min - F_max) / (P_max - P_min)
    print("m ", m)
    # breakpoint()
    # m = 1.5*m
    b = F_max - m*P_max
    # b = F_star - m*Psa_u_star
    F = m*dPsa + b
    
    return F

def get_exp_heart_rate(dPsa, Psa_star, F_star):
    """returns new F

    Args:
        P_sa (_type_): _description_
        P_sa_star (_type_): _description_
        F_prev (_type_): _description_
    """

    d = (F_star - np.exp(-Psa_star))
    # c = np.log(F_star - d) + Psa_star
    F = np.exp(-dPsa) + d

    return F

def get_lower_peripheral_resistance(dPsa, Psa_u_star, Rsl_star, 
                                    Rsl_min, Psa_u_min):
    """ return the target lower peripheral resistance to account for a drop or
    increase in upper systemic blood pressure

    Rsl_min = 100 dynes
    Rsl_max = 200 dynes
    Args:
        dPsa (_type_): _description_
        Psa_u_star (_type_): _description_
        Rsl_star (_type_): _description_
        Rsl_min (_type_): _description_
        Psa_u_min (_type_): _description_
    """
    m = (Rsl_min - Rsl_star) / (Psa_u_star - Psa_u_min)
    b = Rsl_star - m * Psa_u_star
    Rs_l = dPsa * m + b

    return Rs_l