clear all
close all

Psa_u_star = 53*1333; %100mmHg * [ dynes/(mmHg*cm^2) ]
Psa_u = Psa_u_star; 
dP_RA = 2*1333; %2 mm Hg * [ dynes/(mmHg*cm^2) ] %RA

height = 167.64; %cm
Hu = 37.25; %0.5*32;%cm %upper compartment   32
Hl = -(0.5*45);%cm %lower 42

rho = 1; %g/cm^3   density of blood 
g_earth = 980; %gravitational acceleration cm/s^2
G = linspace(0,980*10, 10000); %cm/s^2

Rs = (20)*1333/(1000/60); %default 17.86, **16.49**-20.02 [mmHg *min/L] (too high) systemic resistance (mmHg/(liters/s))
Gs = 1/Rs; %conductance = 1/resistance, resitance = pressure/flow
Gs_u = (1/3)*Gs;
Gs_l = (2/3)*Gs; 
Rs_l = 1/Gs_l;
Rs_u = 1/Gs_u; 
Rp = (1.61*1333)/(1000/60); %pulmonic resistance (mmHg/(liters/s))

C_RVD = (0.0350/1333)*1000; %right-ventricular diastolic compliance (liters/mmHg)
C_LVD = (0.00583/1333)*1000; %left-ventricular diastolic compliance (liters/mmHg)

% 5 risk factors, 2 risk factors, >65, younger
Csa_vec = [0.0012*0.58, 0.0012*0.80, 0.0012, 0.0018, 0.002]

Csa_l = (2/3)*Csa_vec/1333*1000;
Csa_u = (1/3)*Csa_vec/1333*1000;
%Csa_l = [(2/3)*(0.0012/1333)*1000 (2/3)*(0.0018/1333)*1000 (2/3)*(0.002/1333)*1000];
%Csa_u = [(1/3)*(0.0012/1333)*1000 (1/3)*(0.0018/1333)*1000 (1/3)*(0.002/1333)*1000]; 
Csv_l = (2/3)*(0.09/1333)*1000;
Csv_u = (1/3)*(0.09/1333)*1000;
Cs_l = Csa_l + Csv_l;

Cpa = (0.00412/1333)*1000;  %compliance pulm art
Cpv = (0.01/1333)*1000;     %compliance pulm vein
Cp = Cpa + Cpv;
Vtotal = 3.7*1000; %3.7 L * 1000 -> cm^3 


Gs = 1/Rs_u + 1/Rs_l; 
Gs_l = 1/Rs_l;
Ts = Csa_u*Rs_u; %Csa_l*Rs_l
Tp = Rp*Cpa;
Csa = Csa_u+Csa_l;

%P_thorax = linspace(-4*1333,3*1333, 8); %mmHg * dynes/(mmHg*cm^2)
P_thorax = -4.3*1333; %mm Hg -> dynes/cm^2 (P=F/A)
P_RA = P_thorax+dP_RA;

Vd_total_vec = zeros(1,length(G)); 
sol_Vd_Csa_G = zeros(length(P_thorax), length(G)); 
 sol_Q_Csa_G = zeros(length(P_thorax), length(G)); 
 sol_F_Csa_G = zeros(length(P_thorax), length(G)); 
 sol_Ppa_Csa_G = zeros(length(P_thorax),length(G));
 
for j = 1:length(Csa)
    for i = 1:length(G)
        if P_thorax <= -dP_RA %CASE 1
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) - ...
            (Tp*Gs+Csa(j))*Psa_u_star - (Tp*Gs_l+Csa_l(j))*rho*G(i)*Hu...
            - Cs_l(j)*rho*G(i)*(-Hl);
            Q = ((1/Rs_u)+(1/Rs_l))*Psa_u_star+ rho*G(i)*Hu/Rs_l;
            F = Q/(C_RVD*(P_RA-P_thorax));
            Ppv = P_thorax + (C_LVD/C_RVD)*dP_RA;
            Ppa = Ppv+Q*Rp;
            
        elseif P_thorax > -dP_RA && P_thorax < rho*G(i)*Hu - dP_RA %CASE 2
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA - ...
            (Tp*Gs+Csa(j))*Psa_u_star - (Tp*Gs_l+Csa_l(j))*rho*G(i)*Hu ...
            - Cs_l(j)*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs_l)*(P_thorax+dP_RA);
            
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
            Vd_total = Vtotal - Cp*(C_RVD/C_LVD)*dP_RA ...
            - (Tp*Gs+Csa(j))*Psa_u_star - (Tp*Gs+Csa_l(j)-Csv_u)*rho*G(i)*Hu...
            - Cs_l(j)*rho*G(i)*(-Hl) - (Csv_l-Tp*Gs)*(P_thorax+dP_RA);
        
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
        sol_Vd_Csa_G(j,:) = Vd_total_vec; 
        sol_Q_Csa_G(j,:) = Q_vec; 
        sol_F_Csa_G(j,:) = F_vec; 
        sol_Ppa_Csa_G(j,:) = Ppa_vec;
end

%conversions:
G = G/100/(g_earth/100); 
sol_Q_Csa_G = sol_Q_Csa_G*60/1000; 
sol_F_Csa_G = sol_F_Csa_G*60;
sol_Ppa_Csa_G = sol_Ppa_Csa_G/1333;
sol_Vd_Csa_G = sol_Vd_Csa_G/1000;

%%% PLOT OPTIONS %%%
 
%%plot(x, y,myLineColorPref,'color', myLineColorVec,'LineWidth', myLineWidth) %buffer times in black
% myLineColorPref="k-"; %black line
% myLineColorVec = [0,0,0]; %black -> shades of gray
% myLineWidth = 1;
% myLabelFontSize=18;
% alpha=0.1; %increment for gray
% beta=1;

%label_P_thorax = num2str([-3:4]'); %in mm Hg, from -3 to 4

%% saveas(h_overlay_nonDM, myPlotOverlay_nonDM, 'pdf')
%set(gca, 'FontSize', myLabelFontSize)
%print(myfig,figName,"-dpdf", "-S4016,2362")
              

%family of plots for Reserve Volume vs. G for different Pthorax values: 
% 
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
% % xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
% % xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
% 
% h1=figure(101)
% clf(101)
% for i = 1: length(P_thorax)
% 
% plot(G,sol_Q_Pthorax_G(i,:))
% xlabel('G', 'interpreter', 'latex')
% ylabel('Cardiac Output (L/min)', 'interpreter', 'latex')
% hold on
% end
% % xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
% % xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
% 
% 
% h2=figure(102)
% clf(102)
% hold on
%     for i = 1:length(P_thorax)
%         myLineColor=myLineColorVec+alpha*(i-1); %set grayscale
%         plot(G,sol_F_Pthorax_G(i,:), myLineColorPref,'color', myLineColor,'LineWidth', myLineWidthPlot)
%         plot([1,2.2, 3.5],[73,77,100], "k*", 'MarkerSize', 14, 'LineWidth', myLineWidthPlot)
%         plot([3.5],[112], "r*", 'MarkerSize', 14, 'LineWidth', myLineWidthPlot)
% 
%     end
% hold off
% xlabel('G')
% ylabel('Heart Rate (per minute)')
% set(gca, 'FontSize', myLabelFontSize)
% hold on
% 
% % xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
% % xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
% 
% 
% 
h3=figure(103)
clf(103)
% for i = 1:length(P_thorax)
% 
% plot(G,sol_Ppa_Pthorax_G(i,:))
% xlabel('G', 'interpreter', 'latex')
% ylabel('Pulmonary Arterial Pressure (mmHg)', 'interpreter', 'latex')
% hold on
% end
% xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
% xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})

    figure(i)
    hold on
for i = 1:length(Csa)

    plot(G, sol_F_Csa_G)
    xlabel('G')
    ylabel('Heart Rate (per minute)')
    hold on
end

% xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
% % xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})

colorvec={[0,0,0]+0.8, [0,0,0]+0.4, [0,0,0], [0,0,0]+0.4, [0,0,0]+0.8}

figure(100)
for i = 1:length(Csa)
% <<<<<<< HEAD:matlab_controller_2021/varyG_conventional_units.m
%     plot(G, sol_Vd_Csa_G(i,:), 'linewidth', 2, 'Color', [0, 0, 0] + 0.25*i)
% =======
    plot(G, sol_Vd_Csa_G(i,:), 'LineWidth', 2, 'Color', colorvec{i})
%>>>>>>> 16da254 (made plots with varying G and Csa in grayscale for AsMA presentation May 2022):matlab_controller_2021/varyGandCsa_conventional_units.m
    xlabel('G')
    ylabel('Reserve Volume (L)')
    title('Varying Arterial Compliance')
    hold on
end
%<<<<<<< HEAD:matlab_controller_2021/varyG_conventional_units.m
%legend('Csa = 1.2 mL/mmHg','Csa = 1.8 mL/mmHg','Csa = 2.0 mL/mmHg')
%set(gca, 'fontSize', 18)
%=======
legend('58% Csa_{>65}', '80% Csa_{>65}', 'Csa_{>65} = 1.2 mL/mmHg','Csa_{21-29} = 1.8 mL/mmHg','Csa_{50-64} = 2.0 mL/mmHg')
set(gca, 'FontSize', 18)
%>>>>>>> 16da254 (made plots with varying G and Csa in grayscale for AsMA presentation May 2022):matlab_controller_2021/varyGandCsa_conventional_units.m
