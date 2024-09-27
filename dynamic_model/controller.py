import numpy as np

def get_linear_heart_rate(dPsa, F_star, F_min, P_star, P_min):
    """return a linear function of heart rate F in terms of pressure

    Args:
        F_star (_type_): maximum viable pressure for individual
        F_min (_type_): minimum viable pressure 
        P_star (_type_): target / maximum viable upper (coronary) 
            arterial pressure for individual
        P_min (_type_): minimum viable upper (coronary) 
            arterial pressure for individual
        assumes F* > F min
    """
    
    m = (F_min - F_star) / (P_max - P_min)
    b = F_star - m*P_max
    F = m*dPsa + b
    if F > 200/60: 
        F = 200/60 
    if F < 0: 
        F = 0 
       
    
    return F

def get_exp_heart_rate(dPsa, F_star, F_min, Psa_star):
    """returns new F

    Args:
        P_sa (_type_): _description_
        P_sa_star (_type_): _description_
        F_prev (_type_): _description_
    """

    d = (F_star - np.exp(-Psa_star))
    # c = np.log(F_star - d) + Psa_star
    F = np.exp(-dPsa) + d

    # if F > 200/60: 
    #     print("max HR exceeded")
    #     F = 200/60
    # if F < F_min:
    #     print("min viable HR")
    #     F = F_min

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

def get_reserve_venous_volume(dPsa, Psa_u_star, Vsv0_star):
    """_summary_

    Args:
        dPsa_u (_type_): _description_
        dP_RA (_type_): determines distension of RA and Stroke Volume

    Returns
        Vsv0 : lower systemic reserve volume
    """
    m = 1000
    b = Vsv0_star - m*Psa_u_star
    Vsv0 = dPsa*m + b

    print("Vsv0 ", Vsv0, " dPsa ", dPsa, " Psa_u_star ", Psa_u_star,
          "Vsv0_star ", Vsv0_star)
    return Vsv0

