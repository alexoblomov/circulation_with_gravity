%% Model for the steady state circulation with gravity
% written by Alanna Kennard (NYU) and Dr. Charles Peskin (NYU)

% This is a modified verion of the model written by Charles S. Peskin
% adapted from Guyton's model [2], to incorporate the effects of gravity.
% In this model we do not explicity consider the effects of gravity on
% the pulmonary circulation, but do incorporate its effects explicity in 
% the systemic circulation. The lumped systemic compartment is split into
% two compartments for the systemic circulation: above the heart and below 
% the heart. We consider this to be roughly a 1/3 to 2/3 split and assign 
% the arbitrary values of 0.5m as the average distance above the heart and
% 1m as the average distance below the heart.  

% to simulate this in the case that the systemic venous pressure of the 
% upper body is constant and equal to 0, and the systemic venous pressure 
% of the lower body increasing with gravity - which corresponds to a partial 
% collapse of the veins set "simulate" equals to "case 1".

% to simulate this in the case that the systemic venous pressure of the 
% upper body is constant and equal to 0, and the systemic venous pressure 
% of the lower body increasing with gravity and the pressure in the right 
% atrium - which corresponds to a length of partially collapsed veins
% proprotional to the pressure in the right atrium, set "simulate" equals to
% "case 2".

% to simulate this in the case that the systemic venous pressure of the 
% upper body is proportional to the pressure in the right atrium and
% decreasing with gravity, and and the systemic venous pressure 
% of the lower body increasing with gravity and the pressure in the right 
% atrium - which corresponds to a situation with no partial collapse of the
% veins entering the heart, set "simulate" equals to
% "case 3".

% The derivation of the equations for the original model can be found in
% the supplementary materials folder in the file:
% "control_of_the_circulation.pdf"
% The new model can be found in the file:
% "control_of_the_circulation_continued.pdf"
% The equations for the new model can be found in the file:
% "circulation_with_gravity.pdf"

%[1] Guyton A.C.: Circulatory Physiology: Cardiac 
% Output and its Regulation. Saunders, Philadelphia, PA, 1963.

% This program makes use of the labels toolbox written by Chad Greene.
% https://www.mathworks.com/matlabcentral/fileexchange/47421-label

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

%% iterate over the cases and values for Pthorax
cases = {"case 1","case 2","case 3"};
% we iterate over values for pressure in the thorax to find the range of 
% pressures for which each case is true.
Pthorax_range = -10:.02:71; % pressure in the thorax. The flows become negative when Pthorax > 71
solutions = zeros(3,length(Pthorax_range),3); % initialize solutions array.
% this holds the triplet pressure, upper flow and lower flow, for each
% pressure value in the thorax and for each case.

pressures = zeros(3,length(Pthorax_range),4); % this holds the quartet of 
% upper and lower systemic pressures for each value of pressure value in 
% the thorax and for each case

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

for j = 1:length(cases)
simulate = cases{j};
 
for i= 1:length(Pthorax_range)
    Pthorax =  Pthorax_range(i);
% new in model
    g = 980; %cm/s^2
    rho = 1; %g/cm^3
    H_upper = 50; H_lower = -100; %cm

% conversions
    F = F_original / 60; % in seconds
    Pthorax = Pthorax * mmh_to_dynes; % dynes/cm^2 
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
%% wrong row column order in A. uncomment to get plots that are currently in draft
% K = Cpa*Rp+(Cpa+Cpv)/(F*CLVD); % constant term present in each of the three cases
% if (simulate == "case 1")
%     % equations
%     b = [VT-VT0+rho*g*H_lower*(Csa_U+Csv_L);...
%         rho*g*H_upper];
%     A = [K+Csa_U*Rs_U,K+Csa_L*Rs_L;-Rs_U Rs_L]; 
% elseif (simulate == "case 2")
%     % equations
%     b = [VT-VT0-(Csa_U+Csv_L)-(Pthorax-rho*g*H_lower)*(Csv_L+Csa_L);...
%         rho*g*H_upper-Pthorax];
%     A = [K+Csa_U*Rs_U+(Csa_L+Csv_L)/(F*CRVD),K+Csa_L*Rs_L+(Csa_L+Csv_L)/(F*CRVD);...
%         1/(F*CRVD)-Rs_U Rs_L-1/(F*CRVD)];
% elseif (simulate == "case 3")
%     % equations
%       b = [VT-VT0-(Csa_L+Csv_L)*(Pthorax-rho*g*H_lower)-...
%           (Csa_U+Csv_U)*(Pthorax-rho*g*H_upper);0];
%       A = [K+Csa_U*Rs_U+(Csa_U+Csa_L+Csv_U+Csv_L)/(F*CRVD),...
%           K+Csa_L*Rs_L+(Csa_U+Csa_L+Csv_U+Csv_L)/(F*CRVD);-Rs_U Rs_L];
% end
%% correct row column orders. very different plots
K = Cpa*Rp+(Cpa+Cpv)/(F*CLVD); % constant term present in each of the three cases
alpha = (Cpa+Cpv)/CLVD;
dPra = Pra-Pthorax;
if Pthorax <= -dPra %CASE 1
    % equations
    b = [VT-VT0+rho*g*H_lower*(Csa_U+Csv_L);...
        rho*g*H_upper];
    A = [K+Csa_U*Rs_U, -Rs_U; K+Csa_L*Rs_L, Rs_L]; 
elseif Pthorax > (-dPra && Pthorax) < (rho*g*H_upper - dPra) %CASE 2
    % equations
    b = [VT-VT0-(Pthorax-rho*g*H_lower)*(Csv_L+Csa_L);...
        rho*g*H_upper-Pthorax];
    A = [K+Csa_U*Rs_U+(Csa_L+Csv_L)/(F*CRVD), 1/(F*CRVD)-Rs_U ;...
         K+Csa_L*Rs_L+(Csa_L+Csv_L)/(F*CRVD), Rs_L-1/(F*CRVD)];
elseif  Pthorax >= (rho*g*H_upper - dPra)
    % equations
      b = [VT-VT0-(Csa_L+Csv_L)*(Pthorax-rho*g*H_lower)-...
          (Csa_U+Csv_U)*(Pthorax-rho*g*H_upper);0];
      A = [Csa_U*Rs_U+ Cpa*Rp+ alpha/F + (Csa_U+Csa_L+Csv_U+Csv_L)/(F*CRVD),...
          -Rs_U;...
           Csa_L*Rs_L+ alpha/F + (Csa_U+Csa_L+Csv_U+Csv_L)/(F*CRVD), Rs_L];
end
% we solve:
Q = A\b;
%% I. upper and lower flows (Q_U and Q_L)
Q_U = Q(1);
Q_L = Q(2);
%% II. solutions for the other model variables: volumes, pressures in terms of the flows
% in the heart
Pla = Pthorax + (Q_U+Q_L)/CLVD;
Vstroke = CLVD*(Pla- Pthorax);
%Vstroke_star(i) = Vstroke; % get case 2 and 3 values
Pra = Pthorax + (Q_U+Q_L)/CRVD;
if simulate == "case 1"
    Psv_U = 0; 
    Psv_L = -rho*g*H_lower;
elseif((simulate == "case 2" || simulate == "case 3"))
    % get values of systemic venous pressure
    if (simulate == "case 2")
        Psv_U = 0;
        Psv_L = Pra - rho*g*H_lower;
    else
        Psv_U = Pra - rho*g*H_upper;
        Psv_L = Pra - rho*g*H_lower;
    end
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

if (simulate == "case 1")
    if(Pra> Pthorax && Pra < 0)
        solutions(1,i,:) = [Pthorax,Q_U,Q_L];         % index by case, pthoraxvalue, q value
        pressures(1,i,:) = [Psv_U,Psv_L,Psa_U,Psa_L];
        eqn_23 = 1;
        if check
            c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
                Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
                Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
                Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol,simulate); % check if original equations are satisfied by this solution
             consistencies_case1 = consistencies_case1 + c;
            if (A*inv(A)-eye(2) > tol)
                singular_case_1 = singular_case_1+1; 
            end
        end
    else
        solutions(1,i,:) = nan;
        eqn_23 = 0;
    end
elseif (simulate == "case 2")
    if(Pra < rho*g*H_upper && Pra > 0)
        %fprintf("Pra is positive and has value %.2f and is less than rho*g*H_upper = %.2f",Pra,rho*g*H_upper);
        solutions(2,i,:) = [Pthorax,Q_U,Q_L]; % save valid values of Pthorax for the case
        pressures(2,i,:) = [Psv_U,Psv_L,Psa_U,Psa_L];
        eqn_23 = 1;
        if check
            c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
                Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
                Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
                Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol,simulate); % check if original equations are satisfied by this solution
            consistencies_case2 = consistencies_case2 + c;
            s = det(A);
            if (A*inv(A)-eye(2) > tol)
                singular_case_2 = singular_case_2+1;
            end
        end
    else
        % > rho*g*H_upper
        solutions(2,i,:)= nan;
        eqn_23 = 0;
    end
elseif (simulate == "case 3")
    if(Pra > rho*g*H_upper && Pra > 0)
        %fprintf("Pra is positive and has value %.2f and is greater than rho*g*H_upper = %.2f",Pra,rho*g*H_upper);
        solutions(3,i,:) = [Pthorax,Q_U,Q_L]; % save valid values of Pthorax for the case
        pressures(3,i,:) = [Psv_U,Psv_L,Psa_U,Psa_L];
        eqn_23 = 1;
        if (check)
            c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
                Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
                Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
                Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol,simulate); % check if original equations are satisfied by this solution
             consistencies_case3 = consistencies_case3 + c;
            s = det(A);
            if (A*inv(A)-eye(2) > tol)
                singular_case_3 = singular_case_3+1;
            end
        end
    else
        solutions(3,i,:) = nan;
        eqn_23 = 0;
    end
end


end
end
%% Plots
C = linspecer(2);

% find range of solutions for case 1
ksi = isnan(solutions(1,:,1,1));
index_range_1_end = find(ksi,1)-1;
% find range of solutions for case 2
lambda = isnan(solutions(2,:,1,1));
index_range_2_c = find(lambda > 0);
for h=1:length(index_range_2_c)
    if(index_range_2_c(h+1)-index_range_2_c(h) > 1)
        index_range_2_end = index_range_2_c(h+1)-1;
        break;
    end
end
% range of solutions for case 3 are the remaining values
figure
% plot upper and lower flows vs Pthorax for case 1
hold on
p1 = plot(solutions(1,1:index_range_1_end,1,1),solutions(1,1:index_range_1_end,2,1),'LineWidth',2,'Color',C(1,:)); % plot Q_U vs Pthorax
p2 = plot(solutions(1,1:index_range_1_end,1,1),solutions(1,1:index_range_1_end,3,1),'LineWidth',2,'Color',C(2,:)); % plot Q_L vs Pthorax

% plot upper and lower flows vs Pthorax for case 2
hold on
p3 = plot(solutions(2,index_range_1_end+2:index_range_2_end,1,1),solutions(2,index_range_1_end+2:index_range_2_end,2,1),'LineWidth',2,'Color',C(1,:)); % plot Q_U vs Pthorax
p4 = plot(solutions(2,index_range_1_end+2:index_range_2_end,1,1),solutions(2,index_range_1_end+2:index_range_2_end,3,1),'LineWidth',2,'Color',C(2,:)); % plot Q_L vs Pthorax
% plot upper and lower flows vs Pthorax for case 3
hold on
p5 = plot(solutions(3,index_range_2_end+1:end,1,1),solutions(3,index_range_2_end+1:end,2,1),'LineWidth',2,'Color',C(1,:)); % plot Q_U vs Pthorax
p6 = plot(solutions(3,index_range_2_end+1:end,1,1),solutions(3,index_range_2_end+1:end,3,1),'LineWidth',2,'Color',C(2,:)); % plot Q_L vs Pthorax

hold on
title('Upper and Lower flows vs Pressure in the thorax')
xlabel('Pressure (dynes/cm^2)')
ylabel('Flow (cm^3/s)')


% line labels
label(p1,'case 1 : upper flow','location','right',...
    'verticalalignment','bottom','fontweight','bold')
label(p2,'case 1 : lower flow','location','right',...
    'verticalalignment','top','fontweight','bold')

label(p3,'case 2','location','center',...
    'slope','fontweight','bold')
label(p4,'case 2','location','center',...
    'slope','fontweight','bold')

label(p5,'case 3','location','center',...
    'slope','fontweight','bold')
label(p6,'case 3','location','center',...
    'slope','fontweight','bold')
hold off

Z =linspecer(3);
% Q = Q upper + Q lower for cases 1,2,3
Q1 = solutions(1,1:index_range_1_end,2,1)+solutions(1,1:index_range_1_end,3,1);
Q2 = solutions(2,index_range_1_end+1:index_range_2_end,2,1)+solutions(2,index_range_1_end+1:index_range_2_end,3,1);
Q3 = solutions(3,index_range_2_end+1:end,2,1)+solutions(3,index_range_2_end+1:end,3,1);

figure
l1 = plot(solutions(1,1:index_range_1_end,1,1),Q1,'LineWidth',2,'Color',Z(1,:)); % plot Q_U + Q_L vs Pthorax
hold on
l2 = plot(solutions(2,index_range_1_end+1:index_range_2_end,1,1),Q2,'LineWidth',2,'Color',Z(2,:)); % plot Q_U + Q_L vs Pthorax
hold on
l3 = plot(solutions(3,index_range_2_end+1:end,1,1),Q3,'LineWidth',2,'Color',Z(3,:)); % plot Q_U+Q_L vs Pthorax
hold on

% graph labels
title('Flow vs Pressure in the thorax')
xlabel('Pressure (dynes/cm^2)')
ylabel('Flow (cm^3/s)')
hold on
% line labels
label(l1,'case 1','location','right',...
    'verticalalignment','top','fontweight','bold')
label(l2,'case 2','location','center',...
    'slope','fontweight','bold')
label(l3,'case 3','location','center',...
    'slope','fontweight','bold')

hold off

figure
% flows v pressure for case I, II, III
% we plot the pressures as function of the flow. These all implicitly
% depend on Pthorax.
% recall that presssures(case,Pthorax value,pressure value) = [Pthorax,Psv_U,Psv_L,Psa_U,Psa_L];
% case I
hold on
plot(solutions(1,1:index_range_1_end,2,1),pressures(1,1:index_range_1_end,3,1),'LineWidth',2,'Color',C(1,:));  %Psa_U v QU
plot(solutions(1,1:index_range_1_end,2,1),pressures(1,1:index_range_1_end,1,1),'LineWidth',2,'Color',C(1,:));  %Psv_U v QU
plot(solutions(1,1:index_range_1_end,3,1),pressures(1,1:index_range_1_end,4,1),'LineWidth',2,'Color',[.9 0 1]); %Psa_L v QL
plot(solutions(1,1:index_range_1_end,3,1),pressures(1,1:index_range_1_end,2,1),'LineWidth',2,'Color',[.9 0 1]); %Psv_L v QL

% case II
hold on
plot(solutions(2,index_range_1_end+1:index_range_2_end,2,1),pressures(2,index_range_1_end+1:index_range_2_end,3,1),'LineWidth',2,'Color',C(1,:));  %Psa_U v QU
plot(solutions(2,index_range_1_end+1:index_range_2_end,2,1),pressures(2,index_range_1_end+1:index_range_2_end,1,1),'LineWidth',2,'Color',C(1,:));  %Psv_U v QU 
plot(solutions(2,index_range_1_end+1:index_range_2_end,3,1),pressures(2,index_range_1_end+1:index_range_2_end,4,1),'LineWidth',2,'Color',[.9 0 1]); %Psa_L v QL
plot(solutions(2,index_range_1_end+1:index_range_2_end,3,1),pressures(2,index_range_1_end+1:index_range_2_end,2,1),'LineWidth',2,'Color',[.9 0 1]); %Psv_L v QL

% case III
hold on
plot(solutions(3,index_range_2_end+1:end,2,1),pressures(3,index_range_2_end+1:end,3,1),'LineWidth',2,'Color',C(1,:));  %Psa_U v QU
plot(solutions(3,index_range_2_end+1:end,2,1),pressures(3,index_range_2_end+1:end,1,1),'LineWidth',2,'Color',C(1,:));  %Psv_U v QU 
plot(solutions(3,index_range_2_end+1:end,3,1),pressures(3,index_range_2_end+1:end,4,1),'LineWidth',2,'Color',[.9 0 1]); %Psa_L v QL
plot(solutions(3,index_range_2_end+1:end,3,1),pressures(3,index_range_2_end+1:end,2,1),'LineWidth',2,'Color',[.9 0 1]); %Psv_L v QL

hold off

% plot Psv_L with Q_L, Psv_U with Q_U
