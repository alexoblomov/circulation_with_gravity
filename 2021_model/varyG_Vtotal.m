clear all
close all

Psa_u_star = 100*1333; %mmHg * [ dynes/(mmHg*cm^2) ]
Psa_u = Psa_u_star; 
dP_RA = 2*1333;

height = 167.64; %cm
Hu = 0.5*32; %32;%cm %upper   32
Hl = -(0.5*42);%cm %lower 42

rho = 1; %g/cm^3
g_earth = 980; %gravitational acceleration cm/s^2
G = linspace(980,980*10, 10000); %cm/s^2

Rs = (16.49)*1333/(1000/60); %default 17.86, 16.49-20.02 (too high) systemic resistance (mmHg/(liters/s))
Gs = 1/Rs;
Gs_u = (1/3)*Gs;
Gs_l = (2/3)*Gs; 
Rs_l = 1/Gs_l;
Rs_u = 1/Gs_u; 
Rp = (1.61*1333)/(1000/60); %pulmonic resistance (mmHg/(liters/s))

C_RVD = (0.0350/1333)*1000; %right-ventricular diastolic compliance (liters/mmHg)
C_LVD = (0.00583/1333)*1000; %left-ventricular diastolic compliance (liters/mmHg)

Csa_l = (2/3)*(0.0018/1333)*1000;
Csa_u = (1/3)*(0.0018/1333)*1000;
Csv_l = (2/3)*(0.09/1333)*1000;
Csv_u = (1/3)*(0.09/1333)*1000;
Cs_l = Csa_l + Csv_l;

Cpa = (0.00412/1333)*1000; 
Cpv = (0.01/1333)*1000; 
Cp = Cpa + Cpv;
Vtotal = [.9*3.7*1000 .95*3.7*1000 .98*3.7*1000 3.7*1000 1.02*3.7*1000 1.05*3.7*1000]; % cm^3 * 1000


Gs = 1/Rs_u + 1/Rs_l; 
Gs_l = 1/Rs_l;
Ts = Csa_u*Rs_u; %Csa_l*Rs_l
Tp = Rp*Cpa;
Csa = Csa_u+Csa_l;

%P_thorax = linspace(-4*1333,3*1333, 8); %mmHg * dynes/(mmHg*cm^2)
P_thorax = -4.3*1333;
P_RA = P_thorax+dP_RA;

Vd_total_vec = zeros(1,length(G)); 
sol_Vd_Vtotal_G = zeros(length(P_thorax), length(G)); 
 sol_Q_Vtotal_G = zeros(length(P_thorax), length(G)); 
 sol_F_Vtotal_G = zeros(length(P_thorax), length(G)); 
 sol_Ppa_Vtotal_G = zeros(length(P_thorax),length(G));
 
for j = 1:length(Vtotal)
    for i = 1:length(G)
        if P_thorax <= -dP_RA %CASE 1
            Vd_total = Vtotal(j) - Cp*(C_RVD/C_LVD)*(dP_RA) - ...
            (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G(i)*Hu...
            - Cs_l*rho*G(i)*(-Hl);
            Q = ((1/Rs_u)+(1/Rs_l))*Psa_u_star+ rho*G(i)*Hu/Rs_l;
            F = Q/(C_RVD*(P_RA-P_thorax));
            Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
            
        elseif P_thorax > -dP_RA && P_thorax < rho*G(i)*Hu - dP_RA %CASE 2
            Vd_total = Vtotal(j) - Cp*(C_RVD/C_LVD)*dP_RA - ...
            (Tp*Gs+Csa(j))*Psa_u_star - (Tp*Gs_l+Csa_l(j))*rho*G(i)*Hu ...
            - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs_l)*(P_thorax+dP_RA);
            
            Psv_l = -rho*G(i)*Hl + P_thorax + dP_RA;
            Psv_u = 0; 
            Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 

            Qs_u = Psa_u_star/Rs_u;
            Qs_l = (Psa_l-Psv_l)/Rs_l;
            Q = Qs_u + Qs_l;
            F = Q/(C_RVD*(dP_RA));
            Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
            
        elseif P_thorax >= rho*G(i)*Hu - dP_RA %CASE 3
            Vd_total = Vtotal(j) - Cp*(C_RVD/C_LVD)*dP_RA ...
            - (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs+Csa_l-Csv_u)*rho*G(i)*Hu...
            - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs)*(P_thorax+dP_RA);
        
            Psv_l = P_thorax + dP_RA + rho*G(i)*(-Hl);
            Psv_u = P_thorax + dP_RA - rho*G(i)*Hu ; 
            Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 

            Qs_u = (Psa_u_star - Psv_u)/Rs_u;
            Qs_l = (Psa_l-Psv_l)/Rs_l;
            Q = Qs_u + Qs_l;
            F = Q/(C_RVD*(dP_RA));
            Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
        end
        if Vd_total > 0 
        Vd_total_vec(i) = Vd_total; 
      
          Q_vec(i) = Q; 
          F_vec(i) = F; 
          Ppa_vec(i) = Ppa; 
          
        else
          Vd_total_vec(i) = NaN; 
      
          Q_vec(i) = NaN; 
          F_vec(i) = NaN; 
          Ppa_vec(i) = NaN;
        end
    end
        sol_Vd_Vtotal_G(j,:) = Vd_total_vec; 
        sol_Q_Vtotal_G(j,:) = Q_vec; 
        sol_F_Vtotal_G(j,:) = F_vec; 
        sol_Ppa_Vtotal_G(j,:) = Ppa_vec;
end
%conversions:
G = G/100/(g_earth/100); 
sol_Q_Vtotal_G = sol_Q_Vtotal_G*60/1000; 
sol_F_Vtotal_G = sol_F_Vtotal_G*60;
sol_Ppa_Vtotal_G = sol_Ppa_Vtotal_G/1333;
sol_Vd_Vtotal_G = sol_Vd_Vtotal_G/1000;

%%% PLOT OPTIONS %%%
 
%%plot(x, y,myLineColorPref,'color', myLineColorVec,'LineWidth', myLineWidth) %buffer times in black
% myLineColorPref="k-"; %black line
% myLineColorVec = [0,0,0]; %black -> shades of gray
% myLineWidth = 1;
% myLabelFontSize=18;
% alpha=0.1; %increment for gray
% beta=1;

% figure(1)
% for i = 1 : length(Vtotal)
%     plot(G, sol_F_Vtotal_G)
%     xlabel('G')
%     ylabel('Heart Rate (per minute)')
%     hold on
% end
   
figure(1)
for i = 1 : length(Vtotal)
    plot(G, sol_Vd_Vtotal_G(i,:), 'linewidth', 2, 'Color', [0, 0, 0] + 0.15*i)
    xlabel('G')
    ylabel('Reserve Volume (L)')
    title('G Tolerance varying total blood volume')
    hold on
end
legend("-10%", "-5%", "-2%", "3.7 L", "+2%", "+5%")
set(gca, 'fontSize', 18)