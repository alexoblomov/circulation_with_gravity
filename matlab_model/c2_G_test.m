%% Model for the steady state circulation with gravity varied
% This simulation is identical to the prototypical model, except that the
% pressure in the thorax is held constant and value of g is varied

% Cases II and III inconsistent. Pra stays negative, when it should go
% positive. Another variable should change to cope with change in gravity.

%% parameters for the original model
    F_original = 80; % heart rate in beats per minute
% in mmmHs
    Pthorax_original = -4; % pressure in the thorax
% mmHs/(Liters/min)
    Rs_original = 17.86 ; % systemic resistance 
    Rp_original = 1.61; % pulmonary resistance
% in L/mmHs
    CRVD_original = 0.35; % compliance in the right ventricle
    CLVD_original = 0.00583; % compliance in the left ventricle
    Csa_original  = 0.00175; % systemic arterial compliance'
    Csv_original  = 0.00175; % systemic venous compliance Csv_original  = 0.09; tried setting equal to arterial
    Cpa_original  = 0.00412; % pulmonary arterial compliance
    Cpv_original  = 0.01; % pulmonary venous compliance
% in L. 
    Vsa0_original = 0.825; % systemic arterial reserve volume
    Vsv0_original = 3.5; % systemic venous reserve volume
    Vpa0_original = 0.1135; % pulmonary arterial reserve volume
    Vpv0_original = 0.18; % pulmonary venous reserve volume
    VT_original = 5; % total volume
%% conversion rates
    mmh_to_dynes = 0.1*13.6*980;
    L_to_cm3 = 1000;
%% parameters for the new model
% _U denotes above the heart
% _L denotes below the heart
% the flows through the arms are ignored in this model

%% model terms and simulation parameters

% we iterate over different g's to see if the system stays
% physiologically consistent.
% lets get the g's
G = load('g_values.mat');
G_range = G.earthgs; %m/s^2
earth_distances = G.planetd;
% range of values for Pthorax in case 2
    load('PthoraxII.mat');
iterate_multiple = 2;

% holding arrays:
    solutions = zeros(length(case2_Pthorax_values),iterate_multiple*length(G_range),3); % initialize solutions array.
    % this holds the triplet: g, upper flow and lower flow, for each case
    pressures = zeros(length(case2_Pthorax_values),iterate_multiple*length(G_range),4); % this holds the quartet of 
    % pressures for each value of g

% we define the tolerance
tol = 10*exp(-03);

% to check that the solutions are consistent with the original equations
consistencies_case1 = 0; % presume none are consistent
consistencies_case2 = 0;
consistencies_case3 = 0;
% to check that the solutions matrix is well conditioned
singular_case_1 = 0;
singular_case_2 = 0;
singular_case_3 = 0;

% and set check to true to run function that check that solutions are
% consistent
check = true;
simulate = "case 2";

for j = 1:length(case2_Pthorax_values)
 
for i= 1:length(G_range)*iterate_multiple
% new in model
    if mod(i,2) == 1
        g = G_range(mod(i,12))*100; %cm/s^2
    else
        g = (G_range(mod((i-1),12))+G_range(mod((i+1),12)))/2*100;
    end
    rho = 1; %g/cm^3
    H_upper = 50; H_lower = -100; %cm

% conversions
    F = F_original / 60; % in seconds
% resistances in dynes/cm^2 / cm^3/second
    Rs = Rs_original * mmh_to_dynes / (L_to_cm3 / 60);
    Rs_U = 1/(1/Rs*(1/3));
    Rs_L = 1/(1/Rs*(2/3));
    Rp = Rp_original * mmh_to_dynes / (L_to_cm3 / 60);
% compliances in cm^3/  dynes/cm^2
    CRVD = CRVD_original * L_to_cm3 / mmh_to_dynes;
    CLVD = CLVD_original * L_to_cm3 / mmh_to_dynes;
    Csa_U = (1/3) * Csa_original * L_to_cm3 / mmh_to_dynes;
    Csa_L = (2/3) * Csa_original * L_to_cm3 / mmh_to_dynes;
    Csv_U = (1/3) * Csv_original * L_to_cm3 / mmh_to_dynes;
    Csv_L = (2/3) * Csv_original * L_to_cm3 / mmh_to_dynes;
    Cpa = Cpa_original * L_to_cm3 / mmh_to_dynes;
    Cpv = Cpv_original * L_to_cm3 / mmh_to_dynes;
% volumes in cm^3
    VT = VT_original * L_to_cm3;
    Vsa0_U = (1/3)* Vsa0_original * L_to_cm3;
    Vsa0_L = (2/3)* Vsa0_original * L_to_cm3;
    Vsv0_U = (1/3)* Vsv0_original * L_to_cm3;
    Vsv0_L = (2/3)* Vsv0_original * L_to_cm3;
    Vpa0 = Vpa0_original * L_to_cm3;
    Vpv0 = Vpv0_original * L_to_cm3;
    VT0 = Vsa0_U + Vsa0_L + Vsv0_U + Vsv0_L + Vpa0 + Vpv0;
%% Cases
K = Cpa*Rp+(Cpa+Cpv)/(F*CLVD); % constant term present in each of the three cases

    Pthorax = case2_Pthorax_values(j) * mmh_to_dynes; % This is the constant value of the pressure in the thorax we use for case II.
    % This is picked as the median value in the range of values of Pthorax
    % for which case II is satistied
    
    % equations
    b = [VT-VT0-(Csa_U+Csv_L)-(Pthorax-rho*g*H_lower)*(Csv_L+Csa_L);...
        rho*g*H_upper-Pthorax];
    A = [K+Csa_U*Rs_U+(Csa_L+Csv_L)/(F*CRVD),K+Csa_L*Rs_L+(Csa_L+Csv_L)/(F*CRVD);...
        1/(F*CRVD)-Rs_U Rs_L-1/(F*CRVD)];
% we solve:
Q = A\b;
%% I. upper and lower flows (Q_U and Q_L)
Q_U = Q(1);
Q_L = Q(2);
%% II. solutions for the other model variables: volumes, pressures in terms of the flows
% in the heart
Pla = Pthorax + (Q_U+Q_L)/CLVD;
Vstroke = CLVD*(Pla- Pthorax);
Pra = Pthorax + (Q_U+Q_L)/CRVD;
% get values of systemic venous pressure
    if (simulate == "case 2")
        Psv_U = 0;
        Psv_L = Pra - rho*g*H_lower;
    else
        Psv_U = Pra - rho*g*H_upper;
        Psv_L = Pra - rho*g*H_lower;
    end
% systemic
Vsv_U = Vsv0_U + Csv_U*Psv_U;
Vsv_L = Vsv0_L + Csv_L*Psv_L;
Psa_U = Psv_U +Q_U*Rs_U;
Psa_L = Psv_L + Q_L*Rs_L;
Vsa_U = Vsa0_U + Csa_U*Psa_U;
Vsa_L = Vsa0_L + Csa_L*Psa_L;
%pulmonary
Ppv = Pthorax + (Q_U+Q_L)/(F*CLVD);
Ppa = Pthorax + (Q_U+Q_L)*(Rp + 1/(F*CLVD));
Vpv = Vpv0 + Cpv*(Ppv-Pthorax);
Vpa = Vpa0 + Cpa*(Ppa-Pthorax);

%% consistency of solutions
% Here we check that the pressure in the right atrium is within the range
% it should be in order to verify the model is consistent.

    if(Pra < rho*g*H_upper && Pra > 0)
        solutions(j,i,:) = [g,Q_U,Q_L]; % save valid triplet for the case
        pressures(j,i,:) = [Psv_U,Psv_L,Psa_U,Psa_L];
        eqn_23 = 1;
        if check
            c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
                Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
                Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
                Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol,simulate); 
            consistencies_case2 = consistencies_case2 + c;
            s = det(A);
            if (A*inv(A)-eye(2) > tol)
                singular_case_2 = singular_case_2+1; 
            end
        end
    else
        % > rho*g*H_upper
        solutions(j,i,:)= nan;
        eqn_23 = 0;
    end

end
end

%% Plots
% we plot the flows and pressures as function of g.
C = linspecer(iterate_multiple*length(G_range)); %colorscheme

[~,ind] = find(~isnan(solutions(:,1,1,1)));
figure

hold on
title('g vs range of acceptable Pthorax values for case II')
xlabel('Pressure (dynes/cm^2)')
ylabel('g (m/s^)')
hold on
%    solutions = zeros(length(case2_Pthorax_values),iterate_multiple*length(G_range),3); % initialize solutions array.

for i = 1:2*length(G_range)
    ksi = squeeze(isnan(solutions(:,i,1)));
    if sum(ksi) == 0 % all are true
        plot(case2_Pthorax_values,solutions(:,i,1,1),'LineWidth',2,'Color',C(i,:));
        hold on
    elseif sum(ksi) == length(ksi) % none are true
            continue
    else
            range_end = find(ksi,1)-1;
            height = repmat(solutions(:,i,1,1),1,length(case2_Pthorax_values(1:range_end))); % y value
            plot(case2_Pthorax_values(1:range_end),height,'LineWidth',2,'Color',C(i,:));
            hold on
    end
end
hold off
figure

% find range of valid solutions for this case:

% plot upper and lower flows vs g for case 1
hold on
Qu = plot(solutions(1,:,1,1)./100,solutions(1,:,2,1),'LineWidth',2,'Color',C(1,:)); % plot Q_U vs Pthorax
Ql = plot(solutions(1,:,1,1)./100,solutions(1,:,3,1),'LineWidth',2,'Color',C(2,:)); % plot Q_L vs Pthorax

hold on
title('Upper and Lower flows vs g')
xlabel('distance from surface (m)')
ylabel('Flow (cm^3/s)')
set(gca, 'Ticklength', [0 0])
xticklabels(earth_distances)
hold off

figure
% Q = Q upper + Q lower
Q1 = solutions(1,:,2,1)+solutions(1,:,3,1);

plot(solutions(1,:,1,1)./100,Q1,'LineWidth',2,'Color',C(3,:)); % plot Q_U + Q_L vs Pthorax
hold on
title('Flow vs g')
xlabel('distance from surface (m)')
ylabel('Flow (cm^3/s)')
set(gca, 'Ticklength', [0 0])
xticklabels(earth_distances)
hold off

figure
% plot pressures vs g for case 1
hold on
p1 = plot(solutions(1,:,1,1)./100,pressures(1,:,3,1),'LineWidth',2,'Color',C(3,:));  %Psa_U v BV
p2 = plot(solutions(1,:,1,1)./100,pressures(1,:,1,1),'LineWidth',2,'Color',C(1,:));  %Psv_U v BV

p3 = plot(solutions(1,:,1,1)./100,pressures(1,:,4,1),'LineWidth',2,'Color',C(4,:)); %Psa_L v BV
p4 = plot(solutions(1,:,1,1)./100,pressures(1,:,2,1),'LineWidth',2,'Color',C(2,:)); %Psv_L v BV

% line labels
label(p3,'lower systemic arterial','location','right',...
    'verticalalignment','bottom','fontweight','bold')
label(p4,'lower systemic venous','location','right',...
    'verticalalignment','bottom','fontweight','bold')

label(p1,'upper systemic arterial','location','right',...
    'verticalalignment','bottom','fontweight','bold')
label(p2,'upper systemic venous','location','right',...
    'verticalalignment','bottom','fontweight','bold')
% graph labels
title('Pressure vs g')
ylabel('Pressure (dynes/cm^2)')
xlabel('distance from surface (m)')
set(gca, 'Ticklength', [0 0])
xticklabels(earth_distances)
hold off