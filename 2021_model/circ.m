clear all
%this is code to simulate case 1: 
%Pthorax<P_RA<0 
%CGS units

Psa_u_star = 100*1333; %dynes/cm^2
Psa_u = Psa_u_star; 
dP_RA = 2*1333;
height = 68; %cm
Hu = (1/(3*2))*height;%upper  
Hl = -(2/(3*2))*height;%lower

rho = 1; %g/cm^3
G = 980; %acceleration of gravity, Earth (cm/s^2) 

P_thorax = -4*1333; %mmHg * dynes/(mmHg*cm^2)


Rs = (17.86)*1333/(1000/60); %systemic resistance (mmHg/(liters/minute))
Gs = 1/Rs;
Gs_u = (1/3)*Gs;
Gs_l = (2/3)*Gs; 
Rs_l = 1/Gs_l;
Rs_u = 1/Gs_u; 
Rp = (1.61*1333)/(1000/60); %pulmonic resistance (mmHg/(liters/minute))

C_RVD = (0.0350/1333)*1000; %right-ventricular diastolic compliance (liters/mmHg)
C_LVD = (0.00583/1333)*1000; %left-ventricular diastolic compliance (liters/mmHg)

Csa_l = (2/3)*(0.00175/1333)*1000;
Csa_u = (1/3)*(0.00175/1333)*1000; 
Csv_l = (2/3)*(0.09/1333)*1000;
Csv_u = (1/3)*(0.09/1333)*1000;
Cs_l = Csa_l + Csv_l;

Cpa = (0.00412/1333)*1000; 
Cpv = (0.01/1333)*1000; 
Cp = Cpa + Cpv;
Vtotal = 5.0*1000; %cm^3 
 


P_RA = P_thorax+dP_RA;
%Case 1 Assumptions:
Psv_l = -rho*G*Hl;
Psv_u = 0; 

Psa_l = Psa_u_star + rho*G*(Hu-Hl); 

Qs_u = Psa_u_star/Rs_u;
Qs_l = (Psa_l-Psv_l)/Rs_l;
Q = ((1/Rs_u)+(1/Rs_l))*Psa_u_star+ rho*G*Hu/Rs_l;
F = Q/(C_RVD*(P_RA-P_thorax));

Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
Ppa = Ppv+Q*Rp;

Csa = Csa_l+Csa_u;
Gs = 1/Rs_u + 1/Rs_l; 
Gs_l = 1/Rs_l;
Ts = Csa_u*Csa_l;
Tp = Cpa*Rp;
Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - ...
    (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G*Hu - ...
    Cs_l*rho*G*(-Hl);
