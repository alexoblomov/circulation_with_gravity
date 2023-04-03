"""purpose: currently unused. potentially useful for refactor
"""

def ss_controller(P_RA, P_thorax, Vtotal, Cp, C_RVD, C_LVD, Tp, Gs, Csa, 
                  Psa_u_star, Gs_l, Csa_u, Csa_l, Csv_u, Csv_l, rho , g, Hu, Hl,
                  Cs_l, Rs_u, Rs_l, Rp): 
    """controller for steady state model. returns values for:
            - Vtotal (total reserve volume)
            - F (heart rate)
    that keep Psa_upper (upper arterial systemic pressure) and 
              dPRA (pressure differential in the thorax)
    constant (some inital chosen constant values)
    NB: may need to have an epsilon that dpra and psa upper can vary within for
    the time varying model.

    P_RA is P_RA[j] 
    P_thorax is  P_thorax
    """
    dP_RA = P_RA -  P_thorax
    if  P_thorax <= - dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                    - (Tp * Gs + Csa) * Psa_u_star \
                    - (Tp * Gs_l + Csa_l) \
                    * rho * g * Hu - Cs_l * rho * g * (- Hl)
        Q = Gs * Psa_u_star + rho * g * Hu / Rs_l
        F = Q / (C_RVD * (P_RA -  P_thorax))
        return [Q, F, Vd_total]
    
    elif (P_thorax > - dP_RA) and  (P_thorax < rho * g * Hu - dP_RA):
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                    - (Tp * Gs + Csa)* Psa_u_star \
                    - (Tp * Gs_l + Csa_l) * (rho * g * Hu - Cs_l *
                                            rho * g * (- Hl)) \
                    - (Csv_l - Tp * Gs_l) * ( P_thorax + dP_RA)
        Psv_l = - rho * g * Hl +  P_thorax + dP_RA
        Psv_u = 0
        Psa_l = Psa_u_star + rho * g * (Hu - Hl)

        Qs_u = Psa_u_star * 1 / Rs_u # eq 61
        #eq 62
        Qs_l = (Psa_u_star + rho * g * Hu \
                -  P_thorax-dP_RA) * 1/Rs_l
        Q = Qs_u + Qs_l
        F = Q / (C_RVD * (dP_RA))   
        return [Q, F, Vd_total]

    elif  P_thorax >= (rho * g * Hu - dP_RA):
        # print("entered case III")
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - \
                (Tp*Gs + Csa_l + Csa_u)*Psa_u_star \
                +(Tp*Gs + Csa_u - Csv_l) * rho * g * Hu \
                + (Csa_l+Csv_l) * rho * g * (-Hl) \
                + (Csv_l+Csv_u - Tp*Gs)* ( P_thorax - dP_RA)

        Psv_l =  P_thorax + dP_RA + rho * g * (- Hl)
        Psv_u =  P_thorax + dP_RA - rho * g * Hu

        Psa_l = Psa_u_star + rho * g * (Hu - Hl)
        Qs_u = (Psa_u_star - Psv_u) / Rs_u

        Qs_l = (Psa_l - Psv_l) / Rs_l

        Q = Qs_u + Qs_l
        F = Q / (C_RVD * (dP_RA))

        return [Q, F, Vd_total]