clear all
%script to simulate max G value tolerance
Psa_u_star = 100*1333; %dynes/cm^2
Psa_u = Psa_u_star; 
dP_RA = 2*1333;


height = 172; %cm
Hu = (1/(3*2))*height;%upper  
Hl = -(2/(3*2))*height;%lower

rho = 1; %g/cm^3
g_earth = 980; %gravitational acceleration cm/s^2
G = linspace(0,10000, 5000); %cm/s^2

Rs = (17.86)*1333/(1000/60); %systemic resistance (mmHg/(liters/s))
Gs = 1/Rs;
Gs_u = (1/3)*Gs;
Gs_l = (2/3)*Gs; 
Rs_l = 1/Gs_l;
Rs_u = 1/Gs_u; 
Rp = (1.61*1333)/(1000/60); %pulmonic resistance (mmHg/(liters/s))

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
Gs = 1/Rs_u + 1/Rs_l; 
Gs_l = 1/Rs_l;
Ts = Csa_u*Csa_l;
Tp = Rp*Cpa;
Csa = Csa_u+Csa_l;

P_thorax = [-4*1333:1333:4*1333] ; %mmHg * dynes/(mmHg*cm^2)
P_RA = P_thorax+dP_RA;

Vd_total_vec = zeros(1,length(G)); 
sol_Vd_Pthorax_G = zeros(length(P_thorax), length(G)); 
 sol_Q_Pthorax_G = zeros(length(P_thorax), length(G)); 
 sol_F_Pthorax_G = zeros(length(P_thorax), length(G)); 
 sol_Ppa_Pthorax_G = zeros(length(P_thorax),length(G));
 
for j = 1:length(P_thorax)
    for i = 1:length(G)
        if P_thorax(j) <= -dP_RA %CASE 1
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - ...
            (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G(i)*Hu...
            - Cs_l*rho*G(i)*(-Hl);
            Q = ((1/Rs_u)+(1/Rs_l))*Psa_u_star+ rho*G(i)*Hu/Rs_l;
            F = Q/(C_RVD*(P_RA(j)-P_thorax(j)));
            Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
        elseif P_thorax(j) > -dP_RA && P_thorax(j) < rho*G(i)*Hu - dP_RA %CASE 2
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA - ...
            (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G(i)*Hu ...
            - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs_l)*(P_thorax(j)+dP_RA);
            
            Psv_l = -rho*G(i)*Hl + P_thorax(j) + dP_RA;
            Psv_u = 0; 
            Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 

            Qs_u = Psa_u_star/Rs_u;
            Qs_l = (Psa_l-Psv_l)/Rs_l;
            Q = Qs_u + Qs_l;
            F = Q/(C_RVD*(dP_RA));
            Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
        elseif P_thorax(j) >= rho*G(i)*Hu - dP_RA %CASE 3
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA ...
            - (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs+Csa_l-Csv_u)*rho*G(i)*Hu...
            - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs)*(P_thorax(j)+dP_RA);
        
            Psv_l = P_thorax(j) + dP_RA + rho*G(i)*(-Hl);
            Psv_u = P_thorax(j) + dP_RA - rho*G(i)*Hu ; 
            Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 

            Qs_u = (Psa_u_star - Psv_u)/Rs_u;
            Qs_l = (Psa_l-Psv_l)/Rs_l;
            Q = Qs_u + Qs_l;
            F = Q/(C_RVD*(dP_RA));
            Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
        end
        Vd_total_vec(i) = Vd_total; 
          Q_vec(i) = Q; 
          F_vec(i) = F; 
          Ppa_vec(i) = Ppa; 
    end
        sol_Vd_Pthorax_G(j,:) = Vd_total_vec; 
        sol_Q_Pthorax_G(j,:) = Q_vec; 
        sol_F_Pthorax_G(j,:) = F_vec; 
        sol_Ppa_Pthorax_G(j,:) = Ppa_vec;
end

figure(100)
subplot(4,2,1)
plot(G,sol_Vd_Pthorax_G(1,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = -4mmHg', 'Location', 'northeast')
subplot(4,2,2)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(2,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = -3mmHg', 'Location', 'northeast')
subplot(4,2,3)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(3,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = -2mmHg', 'Location', 'northeast')
subplot(4,2,4)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(4,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = -1mmHg', 'Location', 'northeast')
subplot(4,2,5)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(5,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = 0mmHg', 'Location', 'northeast')
subplot(4,2,6)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(6,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = 1mmHg', 'Location', 'northeast')
subplot(4,2,7)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(7,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = 2mmHg', 'Location', 'northeast')
subplot(4,2,8)
title('Reserve Volume varying G')
plot(G,sol_Vd_Pthorax_G(8,:))
xlabel('G, cm/s$^2$', 'interpreter', 'latex')
ylabel('V$^0$, (cm$^3$)', 'interpreter', 'latex')
ylim([0 5000])
legend('Pthorax = 3mmHg', 'Location', 'northeast')

