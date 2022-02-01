function [Q,solution_set_pressures] = original_model_sim(simulate,Pthorax)
% must resolve om for case where Psv = Pra;
%% parameters
    F = 80; % heart rate in beats per minute
% mmHs/(Liters/min)
    Rs = 17.86 ; % systemic resistance 
    Rp = 1.61; % pulmonary resistance
% in L/mmHs
    CRVD = 0.035; % compliance in the right ventricle
    CLVD = 0.00583; % compliance in the left ventricle
    Csa  = 0.00175; % systemic arterial compliance'
    Csv  = 0.09; % systemic venous compliance Csv_original  = 0.09; tried setting equal to arterial
    Cpa  = 0.00412; % pulmonary arterial compliance
    Cpv = 0.01; % pulmonary venous compliance
% in L. 
    Vsa0 = 0.825; % systemic arterial reserve volume
    Vsv0 = 3.5; % systemic venous reserve volume
    Vpa0 = 0.1135; % pulmonary arterial reserve volume
    Vpv0 = 0.18; % pulmonary venous reserve volume
    VT = 5; % total volume
    VT0 = Vsa0 + Vsv0 + Vpa0 + Vpv0;
% the threshold value for the sign of Pra is:
   % T = -((VT-VT0)/CRVD)/(F*(Rs*Csa+Rp*Cpv)+(Cpa+Cpv)/CLVD);
% then if the pressure in the thorax < T, we are in case I or II, if Pthorax > T we are in case II or III    
        
%% Equations

    if ((simulate == "case 1" || simulate == "case 2"))
        Q = (VT-VT0)/(Rs*Csa+Rp*Cpa+(Cpa+Cpv)/(F*CLVD));
        Pra = Pthorax + Q/(F*CRVD);
        Psv = 0;
    elseif (simulate == "case 3" || "case 2")
      %  Psv = Pra;
      Q = ((VT-VT0)-(Csa+Csv)*Pthorax)/...
            (Rs*Csa+Rp*Cpa+(1/F)*((Cpa+Cpv)/CLVD+(Csa+Csv)/CRVD));
        Pra = Pthorax + Q/(F*CRVD);
        Psv = Pra;
    end

Vstroke = Q/F;
Pla = Pthorax + Q/(F*CLVD);
Ppv = Pthorax + Q/(F*CLVD);
Ppa = Pthorax + Q*(Rp+1/(F*CLVD));
Psa = Psv + Q*Rs;

Vpv = Vpv0 + (Cpv*Q)/(F*CLVD);
Vpa = Vpa0 + Q*(Rp*Cpa + Cpa/(F*CLVD));
Vsa = Vsa0 + Csa*Psa;
Vsv = Vsv0 + Csv*Psv;
solution_set_pressures = [Pla,Pra,Ppv,Ppa,Psv,Psa];

end