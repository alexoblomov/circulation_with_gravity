%% Original model controller
%% parameters
simulate = "case 1";
lambda = 0;
% mmHs/(Liters/min)
    Rs = 17.86 ; % systemic resistance 
    Rp = 1.61; % pulmonary resistance
% in L/mmHs
    CRVD = 0.35; % compliance in the right ventricle
    CLVD = 0.00583; % compliance in the left ventricle
    Csa  = 0.00175; % systemic arterial compliance'
    Csv  = 0.00175; % systemic venous compliance Csv_original  = 0.09; tried setting equal to arterial
    Cpa  = 0.00412; % pulmonary arterial compliance
    Cpv  = 0.01; % pulmonary venous compliance
% in L. 
    Vsa0 = 0.825; % systemic arterial reserve volume
    Vsv0 = 3.5; % systemic venous reserve volume
    Vpa0 = 0.1135; % pulmonary arterial reserve volume
    Vpv0 = 0.18; % pulmonary venous reserve volume
    VT = 5; % total volume
% in mmmHs
    if ((simulate == "case 1" || simulate == "case 2"))
        Psv = 0;
        Pthorax = -4;
    elseif (simulate == "case 3" || "case 2")
        Pthorax = 53.82;
    end
% the threshold value for the sign of Pra is:
%    Th = -((VT-VT0)/CRVD)/(F*(Rs*Csa+Rp*Cpv)+(Cpa+Cpv)/CLVD);
% then if the pressure in the thorax < T, we are in case I or II, if Pthorax > T we are in case II or III    

% the starting guess:
Vstroke_star = .07; %L
Q_star = 5.6; % L/m
VT0_star = (0.825+3.5+.1135+.18); %L
F_star = 80; % bpm
Psa_star = 100;  % mmHG
%% Equations
    % constants
    alpha = (Cpa+Cpv)/(CLVD);
    T = Rs*Csa + Rp*Cpa;
    Beta = alpha+(Csa+Csv)/CRVD;
    % the system
    G = @(g0)[((VT-g0(2))-(Csa+Csv)*Pthorax)/(T+Beta/g0(5))-g0(1); ... f1
        lambda*(T*(g0(3)-Psa_star)/Rs+alpha*(g0(4)-Vstroke_star))-g0(2)+VT0_star;... f2 
        Rs*g0(5)*g0(4)-g0(3);... f3 53a
        (VT-g0(2))/(g0(5)*T+alpha)-g0(4);... f4 53b
        ((Psa_star*g0(4))/(g0(3)*Vstroke_star))^lambda-(g1(5)/F_star);... f5
        ];
    % the jacobian of the system
    J = @(g0)[-1,-1/(T+(alpha+(Csa+Csv)/CRVD)/g0(5)),0,0,(Beta*((VT-g0(2))-(Csa+Csv)*Pthorax))/(g0(5)*T+Beta)^2;...df1
        0,-1,lambda*T/Rs,lambda*alpha,0;...df2
        0,0,-1,Rs*g0(5),Rs*g0(4);...df3
        0,-1/(g0(5)*T+alpha),0,-1,(T*(g0(2)-VT))/(g0(5)*T+alpha)^2;... df4
        0,0,lambda*((-Psa_star*g0(4))/(Vstroke_star*g0(3)^2))^(lambda-1),...
        lambda*(Psa_star/(g0(3)*Vstroke_star))^(lambda-1),-1/F_star];
    
    % we choose the initial guess to be close to the steady state values
    a = 0.8; b =1.2;
    epsilon = a+b.*rand(size(G,1));
    
    % the initial guess of solutions
    g0 = [Q_star;VT0_star;...
        Psa_star;Vstroke_star;F_star].*epsilon;
    
    % we solve:
    [S,numIts] = mnewton(G,J,g0,1000);

    if ((simulate == "case 1" || simulate == "case 2"))
        Pra = Pthorax + Q/(F*CRVD);
    elseif (simulate == "case 3" || "case 2")
        Pra = Pthorax + Q/(F*CRVD);
        Psv = Pra;
    end

%% I. solutions for the control variables:

Q       = S(1);% the upper and lower flows (Q_U and Q_L)
VT0     = S(2);% the reserve volume of the systemic veins in the lower body
Psa     = S(3);% the systemic arterial pressure
Vstroke = S(4);% the stroke volume
F       = S(5);% the heart rate

%% solutions for the other model variables
Pla = Pthorax + Q/(F*CLVD);
Ppv = Pthorax + Q/(F*CLVD);
Ppa = Pthorax + Q*(Rp+1/(F*CLVD));

Vpv = Vpv0 + (Cpv*Q)/(F*CLVD);
Vpa = Vpa0 + Q*(Rp*Cpa + Cpa/(F*CLVD));
Vsa = Vsa0 + Csa*Psv+Q*Rs*Csa;
Vsv = Vsv0 + Csv*Psv;

solution_set_volumes = [Vstroke,Vpv,Vpa,Vsa,Vsv];
solution_set_pressures = [Pla,Pra,Ppv,Ppa,Psv,Psa];