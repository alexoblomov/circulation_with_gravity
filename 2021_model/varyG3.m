clear all
%this is code to simulate steady-state case 3: 
%rho*G*Hu < P_RA, no collapse of sys veins
%CGS units

Psa_u_star = 100*1333; %dynes/cm^2
Psa_u = Psa_u_star; 
dP_RA = 2*1333;
P_thorax = 2*1333; %mmHg * dynes/(mmHg*cm^2)
%P_RA = P_thorax+dP_RA
height = 172; %cm
Hu = (1/(3*2))*height;%upper  
Hl = -(2/(3*2))*height;%lower

rho = 1; %g/cm^3
g_earth = 980; %gravitational acceleration cm/s^2
G = linspace(0,10000, 5000); %cm/s^2

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
 
Csa = Csa_l+Csa_u;
Tp = Rp*Cpa;
Csa = Csa_u+Csa_l;
for i = 1 : length(G)

%Case 3:
Psv_l(i) = P_thorax + dP_RA + rho*G(i)*(-Hl);
Psv_u(i) = P_thorax + dP_RA - rho*G(i)*Hu ; 
Psa_l(i) = Psa_u_star + rho*G(i)*(Hu-Hl); 

Qs_u(i) = (Psa_u_star - Psv_u(i))/Rs_u;
Qs_l(i) = (Psa_l(i)-Psv_l(i))/Rs_l;
Q(i) = Qs_u(i) + Qs_l(i);
F(i) = Q(i)/(C_RVD*(dP_RA));

Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
Ppa(i) = Ppv+Q(i)*Rp;

Vd_total(i) = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA ...
    - (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs+Csa_l-Csv_u)*rho*G(i)*Hu...
    - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs)*(P_thorax+dP_RA);
end

figure(7)
plot(G,Q, 'linewidth', 2.5)
title('Cardiac Output varying G')
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('Q, (cm$^3$/s)', 'interpreter', 'latex')

figure(8)
plot(G,Vd_total, 'linewidth', 2.5)
title('Reserve Volume varying G')
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')

figure(9)
plot(G,F, 'linewidth', 2.5)
title('Heartrate varying G')
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('F (s$^{-1}$)', 'interpreter', 'latex')
