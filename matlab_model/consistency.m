function c = consistency(Vsa_U,Vsa_L,Vsv_U,Vsv_L,Vpa,Vpv,Vstroke,Q_U,Q_L,...
            Pla,Pra,Ppv,Ppa,Psa_L,Psv_L,Psa_U,Psv_U,Vsa0_U,Vsa0_L,Vsv0_U,...
            Vsv0_L,Vpa0,Vpv0,VT,F,CLVD,Csa_U,Csa_L,Csv_U,Csv_L,Cpa,Cpv,...
            Rs,Rs_L,Rs_U,CRVD,Pthorax,eqn_23,tol,simulate)
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

switch simulate
    case "case 1"
        % check that the systemic, and pulmonary pressures are positive,
        % but that the pressure in the right atrium is negative. This is
        % necessary so that the resistance of the partially collaped veins
        % remain positive. We cannot have the systemic venous pressure
        % negative, as this would imply a general collapse of the systemic
        % veins.
            if((Psv_U && Psv_L && Pla) >= 0) && Pra < 0
                P = 1;
            else
                P = 0;
                disp("Pressures inconsistent")
            end
            if (Vsv_L && Vsv_U && Vpv && Vpa && Vsa_U && Vsa_L) > 0
                V = 1;
            else
                V = 0;
                disp("Volume inconsistent")
            end
        % for case I, (maybe case II as well)
        % Pra, which should be negative
    case "case 2"
        % In the borderline case Psv = Pra = 0
            if(((Psv_U && Psv_L && Pla) >= 0) && Pra > 0 )|| ((Psv_U + Psv_L) == 0 && Pra == 0)
                P = 1;
            else
                P = 0;
                disp("Pressures inconsistent")
            end
            if (Vsv_L && Vsv_U && Vpv && Vpa && Vsa_U && Vsa_L) > 0
                V = 1;
            else
                V = 0;
                disp("Volume inconsistent")
            end
    case "case 3"
        % in case III Pra & Psv > 0. This is consitent with the fact we are
        % not modeling the partial collapse in this case.
            if(Psv_U && Psv_L && Pla && Pra > 0)
                P = 1;
            else
                P = 0;
                disp("Pressures inconsistent")
            end
            if (Vsv_L && Vsv_U && Vpv && Vpa && Vsa_U && Vsa_L) > 0
                V = 1;
            else
                V = 0;
                disp("Volume inconsistent")
            end
end
% for the case with control, we would expect that the steady state values
% with an induced height difference equal those where the height difference 
% is 0, or equivalently the model without gravity
consistent = [eqn_16a,eqn_16b,eqn_17a,eqn_17b,eqn_18,eqn_19,eqn_20a,...
    eqn_20b,eqn_21a,eqn_21b,eqn_21c,eqn_22,eqn_23,eqn_24,eqn_25,P,V];

if(length(consistent) == sum(consistent))
    c = 1;
else
    c = 0;
    disp("inconsistent");
end