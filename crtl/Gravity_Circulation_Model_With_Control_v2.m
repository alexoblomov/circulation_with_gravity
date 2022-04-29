%% Model for the steady state circulation with gravity
% written by Alanna Kennard (NYU) and Dr. Charles Peskin (NYU)

% This is a modified verion of the model written by Charles S. Peskin
% adapted from Guyton's model [2], to study the reactivity of the optimal 
% controller to the incorporation of the effects of gravity.
% In this model we do not explicity consider the effects of gravity on
% the pulmonary circulation, but do incorporate its effects explicity in 
% the systemic circulation. The lumped systemic compartment is split into
% two compartments for the systemic circulation: above the heart and below 
% the heart. We consider this to be roughly a 1/3 to 2/3 split and assign 
% the arbitrary values of 0.5m as the average distance above the heart and
% 1m as the average distance below the heart.  

% This version incorporates the optimal controller designed for the
% original model, and applies this to the first case: partial collapse of
% the systemic veins entering the heart.

% The derivation of the equations for the original model can be found in
% the supplementary materials folder in the file:
% "control_of_the_circulation.pdf"
% The new model can be found in the file:
% "control_of_the_circulation_continued.pdf"
% The equations for the new model can be found in the file:
% "circulation_with_gravity.pdf"
% and in "model_jacobian.pdf"

%[1] Guyton A.C.: Circulatory Physiology: Cardiac 
% Output and its Regulation. Saunders, Philadelphia, PA, 1963.

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

%% iterate over the values for Pthorax

% we iterate over values for pressure in the thorax to find the range of 
% pressures for which each case is true.
Pthorax_range = -4:.02:71; % pressure in the thorax. The flows become negative when Pthorax > 71

% solutions array. this holds the triplet Pressure, upper flow and lower flow, for each
% pressure value and for each case
%solutions = zeros(3,length(Pthorax_range),3); % initialize solutions array.

% we define the tolerance
tol = 10*exp(-03);

% and the positive constant lambda
lambda = 9e+10; 

% to check that the solutions are consistent with the original equations
consistencies_case1 = 0; % presume none are consistent
consistencies_case2 = 0;
consistencies_case3 = 0;


% and set check to true to run function that check that solutions are
% consistent
check = false;

%for i= 1:length(Pthorax_range)
    Pthorax = Pthorax_range(1); %Pthorax = Pthorax_range(i);
%% control parameters

% uncontrolled values in original model
    Q_U_star = (1/3) * 5.6 * L_to_cm3/60;
    Q_L_star = (2/3) * 5.6 * L_to_cm3/60; 
    Vsv0_L_star = (2/3)* Vsv0_original * L_to_cm3; % in cm^3
    Psa_U_star = 133280; %dynes/cm^2
    Vstroke_star = 0.07*L_to_cm3; % in cm^3
    F_star = F_original / 60; % in seconds

%% model constants
    g = 980; %cm/s^2
    rho = 1; %g/cm^3
    H_upper = 50; H_lower = -100; %cm
% conversions
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
    Csa = Csa_original * L_to_cm3 / mmh_to_dynes;
% volumes in cm^3
    VT = VT_original * L_to_cm3;
    Vsa0_U = (1/3)* Vsa0_original * L_to_cm3;
    Vsa0_L = (2/3)* Vsa0_original * L_to_cm3;
    Vsv0_U = (1/3)* Vsv0_original * L_to_cm3;
    Vpa0 = Vpa0_original * L_to_cm3;
    Vpv0 = Vpv0_original * L_to_cm3;
    VT0 = Vsa0_U + Vsa0_L + Vsv0_U + Vsv0_L_star + Vpa0 + Vpv0; % cm^3
    Vsv_L_star = (2/3)*3.5*L_to_cm3;


%% In Case I
% for convenience, we write the lumped parameters
alpha = (Cpa+Cpv)/(F_star*CLVD);
T = Rs_U*Csa_U + Rp*Cpa; T_star =T; Rs_star = Rs;
K = Cpa*Rp+alpha; 

    % the set of equations that describe the system
        G = @(g0)[Vsa0_U+Vsa0_L+Vsv0_U+g0(3)+Vpa0+Vpv0-VT-rho*g*(H_lower)*(Csv_L+Csa_L)+...
        g0(1)*(K)+g0(2)*(Csa_L*Rs_L+K);... %f1 eqn for VT
        g0(2)*Rs_L-g0(1)*Rs_U-rho*g.*H_upper;... %f2 eqn for sum of flows
        lambda*((T_star/Rs_star)*(g0(4)-Psa_U_star)+alpha*(g0(5)-Vstroke_star))-g0(3)+Vsv0_L_star;... %f3 eqn for VsvL reserve volume
        Rs_U*g0(6)*g0(5)-g0(4);... %f4 eqn for Psa_U
        (VT-VT0)/(g0(6)*(Rs*Csa+Rp*Cpa)+alpha)-g0(5);... %f5 eqn for Vstroke
        ((Psa_U_star*g0(5))/(g0(4)*Vstroke_star))^lambda-g0(6)/F_star];%f6 eqn for F
    
    % the jacobian of the system
    J = @(g0)[Rp*Cpa+alpha,Csa_L*Rs_L+K,1,Csa_U,0,0;... % df1
        -Rs_U,Rs_L,0,0,0,0;... %df2 
        0,0,-1,(-lambda*T_star)/Rs_star,lambda*alpha,0;... % df3
        0,0,0,-1,Rs_U*g0(6),Rs_U*g0(5);... %df4
        0,0,-1,lambda*((Vsv_L_star-Vsv0_L_star)*T*F_star)/(Psa_U_star*(T*F_star+alpha)),...
        lambda*(alpha*(Vsv_L_star-Vsv0_L_star))/(Vstroke_star*(T*F_star+alpha)),0;... %df5
        0,0,0,(-lambda/g0(4))*((Psa_U_star*g0(5))/(g0(4)*Vstroke_star))^lambda,...
        (lambda/g0(5))*((Psa_U_star*g0(5))/(g0(4)*Vstroke_star))^lambda,-1/F_star;...% df6
        ]; %df6

    % we choose the initial guess to be close to the steady state values
    a = 0.8; b =1.2;
    epsilon = a+b.*rand(size(G,1));
    
    % then the initial guess for the solution is:
    g0 = [Q_U_star;Q_L_star;Vsv0_L_star;...
        Psa_U_star;Vstroke_star;F_star].*epsilon;
    
    % we solve:
    [S,numIts] = mnewton(G,J,g0,1000);

%% I. solutions for the control variables:

Q_U     = S(1);% the upper and lower flows (Q_U and Q_L)
Q_L     = S(2);
Vsv0_L  = S(3);% the reserve volume of the systemic veins in the lower body
Psa_U   = S(4);% the systemic venous pressure in the upper body
Vstroke = S(5);% the stroke volume
F       = S(6);% the heart rate
 
%% II. solutions for the other model variables: volumes, pressures in terms of the unkowns

% systemic
Vsv_U = Vsv0_U + Csv_U*Psv_U;
Vsv_L = Vsv0_L + Csv_L*Psv_L;
Psa_L = Psv_L + Q_L*Rs_L;
Vsa_U = Vsa0_U + Csa_U*Psa_U;
Vsa_L = Vsa0_L + Csa_L*Psa_L;
%pulmonary
Ppv = Pthorax + (Q_U+Q_L)/(F*CLVD);
Ppa = Pthorax + (Q_U+Q_L)*(Rp + 1/(F*CLVD));
Vpv = Vpv0 + Cpv*(Ppv-Pthorax);
Vpa = Vpa0 + Cpa*(Ppa-Pthorax);

%end
%% consistency of solutions
% Here we check that the pressure in the right atrium is within the range
% it should be in order to verify the model is consistent.
if false

    if(Pra> Pthorax && Pra < 0)
        %fprintf("Pra is negative and has value %.2f and is less than Pthorax = %.2f",Pra,Pthorax);
        solutions(1,i,:) = [Pthorax,Q_U,Q_L];         % index by case, pthoraxvalue, q value
        eqn_23 = 1;
        if check
            c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
                Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
                Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
                Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol); % check if original equations are satisfied by this solution
             consistencies_case1 = consistencies_case1 + c;
        end
    else
        solutions(1,i,:) = nan;
        eqn_23 = 0;
    end
%end


%figure
% find range of solutions for case 1
ksi = isnan(solutions(1,:,1,1));
index_range_1_end = find(ksi,1)-1;
% find range of solutions for case 2
beta = isnan(solutions(2,:,1,1));
index_range_2_c = find(beta > 0);
for h=1:length(index_range_2_c)
    if(index_range_2_c(h+1)-index_range_2_c(h) > 1)
        index_range_2_end = index_range_2_c(h+1)-1;
        break;
    end
end
% range of solutions for case 3 are the remaining values

% plot upper and lower flows vs Pthorax for case 1
hold on
plot(solutions(1,1:index_range_1_end,1,1),solutions(1,1:index_range_1_end,2,1),'LineWidth',2,'Color',[1 0 1]); % plot Q_U vs Pthorax
plot(solutions(1,1:index_range_1_end,1,1),solutions(1,1:index_range_1_end,3,1),'LineWidth',2,'Color',[0 0 1]); % plot Q_L vs Pthorax

% plot upper and lower flows vs Pthorax for case 2
hold on
plot(solutions(2,index_range_1_end+1:index_range_2_end,1,1),solutions(2,index_range_1_end+1:index_range_2_end,2,1),'LineWidth',2,'Color',[1 0 1]); % plot Q_U vs Pthorax
plot(solutions(2,index_range_1_end+1:index_range_2_end,1,1),solutions(2,index_range_1_end+1:index_range_2_end,3,1),'LineWidth',2,'Color',[0 0 1]); % plot Q_L vs Pthorax
% plot upper and lower flows vs Pthorax for case 3
hold on
plot(solutions(3,index_range_2_end+1:end,1,1),solutions(3,index_range_2_end+1:end,2,1),'LineWidth',2,'Color',[1 0 1]); % plot Q_U vs Pthorax
plot(solutions(3,index_range_2_end+1:end,1,1),solutions(3,index_range_2_end+1:end,3,1),'LineWidth',2,'Color',[0 0 1]); % plot Q_L vs Pthorax

hold on
title('Upper and Lower flows vs Pressure in the thorax')
xlabel('Pressure (dynes/cm^2)')
ylabel('Flow (cm^3/s)')
hold off

figure
% Q = Q upper + Q lower for cases 1,2,3
Q1 = solutions(1,1:index_range_1_end,2,1)+solutions(1,1:index_range_1_end,3,1);
Q2 = solutions(2,index_range_1_end+1:index_range_2_end,2,1)+solutions(2,index_range_1_end+1:index_range_2_end,3,1);
Q3 = solutions(3,index_range_2_end+1:end,2,1)+solutions(3,index_range_2_end+1:end,3,1);

plot(solutions(1,1:index_range_1_end,1,1),Q1,'LineWidth',2,'Color',[1 0 1]); % plot Q_U + Q_L vs Pthorax
hold on
plot(solutions(2,index_range_1_end+1:index_range_2_end,1,1),Q2,'LineWidth',2,'Color',[1 0 1]); % plot Q_U + Q_L vs Pthorax
hold on
plot(solutions(3,index_range_2_end+1:end,1,1),Q3,'LineWidth',2,'Color',[1 0 1]); % plot Q_U+Q_L vs Pthorax
hold on

title('Flow vs Pressure in the thorax')
xlabel('Pressure (dynes/cm^2)')
ylabel('Flow (cm^3/s)')
hold off

end
function c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
            Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
            Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
            Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol)
% Here we check that the solutions to these equations solve the original
% set of equations. of Supplementary materials/control of the circulation
% p 17
eqn_16a = (Vsa_U - (Vsa0_U + Csa_U*Psa_U)) < tol;
eqn_16b = (Vsa_L - (Vsa0_L + Csa_L*Psa_L)) < tol;
eqn_17a = (Vsv_U - (Vsv0_U + Csv_U*Psv_U)) < tol;
eqn_17b = (Vsv_L - (Vsv0_L + Csv_L*Psv_L)) < tol;
eqn_18  = (Vpa - (Vpa0 + Cpa*(Ppa - Pthorax))) < tol;
eqn_19  = (Vpv - (Vpv0 + Cpv*(Ppv - Pthorax))) < tol;
eqn_20a = (Vstroke - CLVD*(Pla - Pthorax)) < tol;
eqn_20b = (Vstroke - CRVD*(Pra - Pthorax)) < tol;
eqn_21a = (Q_U+Q_L - F*Vstroke) < tol; 
eqn_21b = (Q_U - (Psa_U-Psv_U)/Rs_U) < tol;
eqn_21c = (Q_L- (Psa_L-Psv_L)/Rs_L) < tol;
eqn_22  = (Ppv - Pla < tol); 
eqn_24  = VT - (Vsa_U + Vsa_L + Vsv_U + Vsv_L + Vpa + Vpv) < tol;
eqn_25  = (1/Rs - (1/Rs_U + 1/Rs_L)) < tol;

consistent = [eqn_16a,eqn_16b,eqn_17a,eqn_17b,eqn_18,eqn_19,eqn_20a,...
    eqn_20b,eqn_21a,eqn_21b,eqn_21c,eqn_22,eqn_23,eqn_24,eqn_25];

if(length(consistent) == sum(consistent))
    c = 1;
else
    c = 0;
    disp("inconsistent");
end
end