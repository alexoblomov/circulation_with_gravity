clear all
close all

Psa_u_star = 100*1333; %mmHg * [ dynes/(mmHg*cm^2) ]
Psa_u = Psa_u_star; 
dP_RA = 2*1333;

height = 167.64; %cm
Hu = 0.5*[30:1:40]; %32%32;%cm %upper   32
Hl = -(0.5*[30:1:50]);%cm %lower 42

rho = 1; %g/cm^3
g_earth = 980; %gravitational acceleration cm/s^2
%G = linspace(980,980*10, 10000); %cm/s^2

Rs = (16.49)*1333/(1000/60); %default 17.86, 16.49-20.02 (too high) systemic resistance (mmHg/(liters/s))
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
Vtotal = 3.7*1000; %3.7 cm^3

Csa = Csa_l+Csa_u;
Gs = 1/Rs_u + 1/Rs_l; 
Gs_l = 1/Rs_l;
Ts = Csa_u*Rs_u; %Csa_l*Rs_l
Tp = Rp*Cpa;
Csa = Csa_u+Csa_l;

G= NaN(length(Hu), length(Hl))

%P_thorax = linspace(-4*1333,-2*1333, 2); %mmHg * dynes/(mmHg*cm^2) %went top 3*13333
P_thorax = [-4]*1333;
P_RA = P_thorax+dP_RA;


%%%%
figure(1)
hold on
for i=1:length(Hl)
    G(:,i) = (Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - (Tp*Gs+Csa)*Psa_u_star)/...
            (Tp*Gs_l+Csa_l)*rho*Hu+ Cs_l*rho*(-Hl(i));
    plot(Hu*2, G(:,i)/g_earth, 'k*-')
    title(num2str(-Hl(i)*2))
end
% xlabel('Upper Height (cm)_')
% ylabel('Maximum G Tolerance')
% hold off

G_new = G/g_earth; 

figure(2)
heatmap(G)
h.XDisplayLabels = num2str((Hu*2)')
h.YDisplayLabels = num2str((Hl*2)')

heatmap(-Hl*2, Hu*2,G_new)
Xlabel("Height Lower (cm)")
Ylabel("Height Upper (cm)")

%%%%%%%%%%%%
G = zeros(length(P_thorax),1);

Vd_total_vec = zeros(1,length(G)); 
sol_Vd_Pthorax_G = zeroxs(length(P_thorax), length(G)); 
 sol_Q_Pthorax_G = zeros(length(P_thorax), length(G)); 
 sol_F_Pthorax_G = zeros(length(P_thorax), length(G)); 
 sol_Ppa_Pthorax_G = zeros(length(P_thorax),length(G));
 
 
 %% Calculate G-tolerace
 %set Vd_total = 0 and solve for G

 

for j = 1:length(P_thorax)
    
    if P_thorax(j) <= -dP_RA 
        G(j) = (Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - (Tp*Gs+Csa)*Psa_u_star)/...
            (Tp*Gs_l+Csa_l)*rho*Hu+ Cs_l*rho*(-Hl);
    
    else
         %CASE 2
        G_case2 = (Vtotal - Cp*(C_RVD/C_LVD)*dP_RA - (Tp*Gs+Csa)*Psa_u_star - (Csv_l-Tp*Gs_l)*(P_thorax(j)+dP_RA))/ ...
         ((Tp*Gs_l+Csa_l)*rho*Hu +Cs_l*rho*(-Hl));

        %Case 3
        G_case3 = (Vtotal - Cp*(C_RVD/C_LVD)*dP_RA - (Tp*Gs+Csa)*Psa_u_star - (Csv_l-Tp*Gs)*(P_thorax(j)+dP_RA))/...
            ((Tp*Gs+Csa_l-Csv_u)*rho*Hu + Cs_l*rho*(-Hl));

        %only save if it meets the requirementys
        if P_thorax(j) > -dP_RA && P_thorax(j) < (rho*G_case2*Hu - dP_RA)
            G(j) = G_case2
        elseif P_thorax(j) >= rho*G_case3*Hu - dP_RA %CASE 3
            G(j) = G_case3
        else
            G(j) = NaN
        end
 
    end 
end

 %% CALCULATE VARIOUS  
% for j = 1:length(P_thorax)
%     for i = 1:length(G)
%         if P_thorax(j) <= -dP_RA %CASE 1
%             Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - ...
%             (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G(i)*Hu...
%             - Cs_l*rho*G(i)*(-Hl);
%             Q = ((1/Rs_u)+(1/Rs_l))*Psa_u_star+ rho*G(i)*Hu/Rs_l;
%             F = Q/(C_RVD*(P_RA(j)-P_thorax(j)));
%             Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
%             Ppa = Ppv+Q*Rp;
%             
%         elseif P_thorax(j) > -dP_RA && P_thorax(j) < rho*G(i)*Hu - dP_RA %CASE 2
%             Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA - ...
%             (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs_l+Csa_l)*rho*G(i)*Hu ...
%             - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs_l)*(P_thorax(j)+dP_RA);
%             
%             Psv_l = -rho*G(i)*Hl + P_thorax(j) + dP_RA;
%             Psv_u = 0; 
%             Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 
% 
%             Qs_u = Psa_u_star/Rs_u;
%             Qs_l = (Psa_l-Psv_l)/Rs_l;
%             Q = Qs_u + Qs_l;
%             F = Q/(C_RVD*(dP_RA));
%             Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
%             Ppa = Ppv+Q*Rp;
%             
%         elseif P_thorax(j) >= rho*G(i)*Hu - dP_RA %CASE 3
%             Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA ...
%             - (Tp*Gs+Csa)*Psa_u_star - (Tp*Gs+Csa_l-Csv_u)*rho*G(i)*Hu...
%             - Cs_l*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs)*(P_thorax(j)+dP_RA);
%         
%             Psv_l = P_thorax(j) + dP_RA + rho*G(i)*(-Hl);
%             Psv_u = P_thorax(j) + dP_RA - rho*G(i)*Hu ; 
%             Psa_l = Psa_u_star + rho*G(i)*(Hu-Hl); 
% 
%             Qs_u = (Psa_u_star - Psv_u)/Rs_u;
%             Qs_l = (Psa_l-Psv_l)/Rs_l;
%             Q = Qs_u + Qs_l;
%             F = Q/(C_RVD*(dP_RA));
%             Ppv = P_thorax(j) + (C_LVD/C_RVD)*dP_RA;
%             Ppa = Ppv+Q*Rp;
%         end
%         if Vd_total > 0 
%         Vd_total_vec(i) = Vd_total; 
%       
%           Q_vec(i) = Q; 
%           F_vec(i) = F; 
%           Ppa_vec(i) = Ppa; 
%           
%         else
%           Vd_total_vec(i) = NaN; 
%       
%           Q_vec(i) = NaN; 
%           F_vec(i) = NaN; 
%           Ppa_vec(i) = NaN;
%         end
%     end
%         sol_Vd_Pthorax_G(j,:) = Vd_total_vec; 
%         sol_Q_Pthorax_G(j,:) = Q_vec; 
%         sol_F_Pthorax_G(j,:) = F_vec; 
%         sol_Ppa_Pthorax_G(j,:) = Ppa_vec;
% end
% 
 %conversions:
 G = G/100/(g_earth/100); 
% sol_Q_Pthorax_G = sol_Q_Pthorax_G*60/1000; 
% sol_F_Pthorax_G = sol_F_Pthorax_G*60;
% sol_Ppa_Pthorax_G = sol_Ppa_Pthorax_G/1333;
%sol_Vd_Pthorax_G = sol_Vd_Pthorax_G/1000;
% 
% %%% PLOT OPTIONS %%%
%  
% %%plot(x, y,myLineColorPref,'color', myLineColorVec,'LineWidth', myLineWidth) %buffer times in black
% myLineColorPref="k-"; %black line
% myLineColorVec = [0,0,0]; %black -> shades of gray
% myLineWidth = 1;
% myLabelFontSize=18;
% alpha=0.1; %increment for gray
% beta=1;
% 
% %label_P_thorax = num2str([-3:4]'); %in mm Hg, from -3 to 4
% 
% %% saveas(h_overlay_nonDM, myPlotOverlay_nonDM, 'pdf')
% %set(gca, 'FontSize', myLabelFontSize)
% %print(myfig,figName,"-dpdf", "-S4016,2362")
%               
% 
% %family of plots for Reserve Volume vs. G for different Pthorax values: 
% % 
% h=figure(100)
% clf(100)
% 
% hold on
% 
%     for i = 1:length(P_thorax)
%         myLineColor=myLineColorVec+alpha*(i-1); %set grayscale
%         %myLineWidthPlot=myLineWidth+beta*(i-1); %set linewidth
%         myLineWidthPlot=myLineWidth; %do not vary line width
%         plot(G,sol_Vd_Pthorax_G(i,:),myLineColorPref,'color', myLineColor,'LineWidth', myLineWidthPlot)
%     end
% 
%     xlabel('G')
%     ylabel('Reserve Volume (L)')
%     set(gca, 'FontSize', myLabelFontSize)
%     
% hold off
% 
% saveas(h, "QvsV0", 'pdf')
% saveas(h, "QvsV0", 'png')
% 
% %% HEART RATE PLOT %%
% h2=figure(101)
% clf(101)
% hold on
%     for i = 1:length(P_thorax)
%         myLineColor=myLineColorVec+alpha*(i-1); %set grayscale
%         plot(G,sol_F_Pthorax_G(i,:), myLineColorPref,'color', myLineColor,'LineWidth', myLineWidthPlot)
%     end
% hold off
% xlabel('G')
% ylabel('Heart Rate (per minute)')
% set(gca, 'FontSize', myLabelFontSize)
% saveas(h, "QvsV0", 'pdf')
% saveas(h, "QvsV0", 'png')
% 
% 
% %% CARDIAC OUTPUT PLOT %%
% h3=figure(102)
% clf(102)
% for i = 1: length(P_thorax)
% 
% plot(G,sol_Q_Pthorax_G(i,:))
% xlabel('G', 'interpreter', 'latex')
% ylabel('Cardiac Output (L/min)', 'interpreter', 'latex')
% hold on
% end

