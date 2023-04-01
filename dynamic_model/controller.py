import numpy as np

def get_heart_rate(dPsa, F_max, F_min, P_max, P_min):
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
    b = F_max - m*P_max
    # b = F_star - m*Psa_u_star
    F = m*dPsa + b
    
    return F

def get_HR(Psa, Psa_star, F_star, F_max, F_min):
    """returns new F

    Args:
        P_sa (_type_): _description_
        P_sa_star (_type_): _description_
        F_prev (_type_): _description_
    """
    # if (P_sa < P_sa_star) and (F < F_max):
    #     F = F_prev + 10
    # elif (P_sa > P_sa_star) and:
    d = (F_star - np.exp(-Psa_star))
    # c = np.log(F_star - d) + Psa_star
    F = np.exp(-Psa) + d

    return F