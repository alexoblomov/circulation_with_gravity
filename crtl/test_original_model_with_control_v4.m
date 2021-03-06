%% Original model controller
% without PSA
%% parameters
simulate = "case 1";
lambda = exp(100); %lambda large
% in mmHs/(Liters/min)
    Rs = 17.86 ; % systemic resistance 
    Rp = 1.61; % pulmonary resistance
% in L/mmHs
    CRVD = 0.035; % compliance in the right ventricle
    CLVD = 0.00583; % compliance in the left ventricle
    Csa  = 0.00175; % systemic arterial compliance
    Csv  = 0.00175; % systemic venous compliance
    Cpa  = 0.00412; % pulmonary arterial compliance
    Cpv  = 0.01; % pulmonary venous compliance
% in L. 
    Vsa0 = 0.825; % systemic arterial reserve volume
    Vsv0 = 3.5; % systemic venous reserve volume
    Vpa0 = 0.1135; % pulmonary arterial reserve volume
    Vpv0 = 0.18; % pulmonary venous reserve volume
    VT = 5; % total volume
% in mmmHs
    if (simulate == "case 1")
        Psv = 0;
        Pthorax = -4;
    elseif (simulate == "case 2")
        Pthorax = 20;
    end
% the starting guess:
Vstroke_star = .07;
Q_star = 5.6;
VT0_star = (0.825+3.5+.1135+.18);
F_star = 80;
%% Equations
    % constants
    alpha = (Cpa+Cpv)/(CLVD);
    T = Rs*Csa + Rp*Cpa;
    Beta = alpha+(Csa+Csv)/CRVD;
    Rs_star = Rs;
    % the system
    G = @(g0)[lambda*(T*(Rs*g0(3)*g0(2)-Rs_star*F_star*Vstroke_star)/Rs+alpha*(g0(2)-Vstroke_star))-g0(1)+VT0_star;... f2 - equation 91
        (VT-g0(1))/(g0(3)*T+alpha)-g0(3);... f4 - equation 53.b
        ((Rs_star*F_star)/(Rs*g0(3)))^lambda-(g0(3)/F_star)];% f5 - equation 85
    % the jacobian of the system
    % dfVT0, dfVstroke, dfF
    J = @(g0)[-1,lambda*T*Rs*g0(3)/Rs+lambda*alpha,lambda*T*Rs*g0(2)/Rs_star;... -d91
        -1/(g0(3)*T+alpha),-1,(T*(g0(1)-VT))/(g0(3)*T+alpha)^2;... -d53.b
        0,0,-lambda/g0(3)*((Rs*F_star)/(Rs*g0(3)))^lambda-1/F_star]; %d85
    
    % we choose the initial guess to be close to the steady state values
    a = 0.8; b =1.2;
    epsilon = a+b.*rand(size(G,1));
    
    % the initial guess of solutions
    g0 = [VT0_star;Vstroke_star;F_star].*epsilon;
    % we solve:
    [S,numIts] = mnewton(G,J,g0,1000);

    if (simulate == "case 1")
        Pra = Pthorax + Q/(F*CRVD);
    elseif (simulate == "case 2")
        Pra = Pthorax + Q/(F*CRVD);
        Psv = Pra;
    end

%% I. solutions for the control variables:

VT0     = S(1);% the reserve volume of the systemic veins in the lower body
Vstroke = S(2);% the stroke volume
F       = S(3);% the heart rate
%% solutions for the other model variables

    if (simulate == "case 1")
        Q = (VT-VT0)/(Rs*Csa+Rp*Cpa+(Cpa+Cpv)/(F*CLVD));
        Pra = Pthorax + Q/(F*CRVD);
        Psv = 0;
    elseif (simulate == "case 2")
        Q = ((VT-VT0)-(Csa+Csv)*Pthorax)/...
            (Rs*Csa+Rp*Cpa+(1/F)*((Cpa+Cpv)/CLVD+(Csa+Csv)/CRVD));
        Pra = Pthorax + Q/(F*CRVD);
        Psv = Pra;
    end
    
Pla = Pthorax + Q/(F*CLVD);
Ppv = Pthorax + Q/(F*CLVD);
Ppa = Pthorax + Q*(Rp+1/(F*CLVD));
Psa = Rs*F*Vstroke;

Vpv = Vpv0 + (Cpv*Q)/(F*CLVD);
Vpa = Vpa0 + Q*(Rp*Cpa + Cpa/(F*CLVD));
Vsa = Vsa0 + Csa*Psv+Q*Rs*Csa;
Vsv = Vsv0 + Csv*Psv;
