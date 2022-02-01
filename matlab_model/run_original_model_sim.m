%% Simulation of original model
%% iterate over the cases and values for Pthorax
cases = {"case 1","case 2"};
% we iterate over values for pressure in the thorax to find the range of 
% pressures for which each case is true.
Pthorax_range = -30:.02:71; % pressure in the thorax. The flows become negative when Pthorax > 71
solutions = zeros(2,length(Pthorax_range),2); % initialize solutions array.
% this holds the triplet pressure, upper flow and lower flow, for each
% pressure value in the thorax and for each case.

% we define the conversion factors
    L_min_to_cm3_s = 1000/60; % liters/min to cm3 per s
    mmh_to_dynes = 0.1*13.6*980;
    L_to_cm3 = 1000;
    
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
 
    for i = 1:length(Pthorax_range)
        Pthorax =  Pthorax_range(i);
    % simulation
    [Q,solution_set_pressures] = original_model_sim(simulate,Pthorax);
    Psv = solution_set_pressures(5);
    Pra = solution_set_pressures(2);
    % conversions
    Q = Q * L_min_to_cm3_s;
    Pthorax = Pthorax * mmh_to_dynes;
    Pra = Pra * mmh_to_dynes;
    Psv = Psv * mmh_to_dynes;
    % consistency of solutions
    % Here we check that the pressure in the right atrium is within the range
    % it should be in order to verify the model is consistent.

    if (simulate == "case 1")
        if(Pra <= 0)
            solutions(j,i,:) = [Pthorax,Q];
        else
            solutions(j,i,:) = nan;
        end
        
    elseif (simulate == "case 2")
        if(Pra >= 0)
            solutions(j,i,:) = [Pthorax,Q];
        else
            solutions(j,i,:) = nan;
        end
    end
    end
end

%% Plots
C = linspecer(4);
figure
% find range of solutions for case 1
ksi = isnan(solutions(1,:,1));
index_range_1_end = find(ksi,1)-1;
% plot flow vs Pthorax for case 1
hold on
p1 = plot(solutions(1,1:index_range_1_end,1),solutions(1,1:index_range_1_end,2),'LineWidth',2,'Color',C(3,:)); % plot Q vs Pthorax
% plot flow vs Pthorax for case 3
% plot upper and lower flows vs Pthorax for case 3
hold on
p2 = plot(solutions(2,index_range_1_end+1:end,1),solutions(2,index_range_1_end+1:end,2),'LineWidth',2,'Color',C(4,:)); % plot Q_U vs Pthorax
hold on

% graph labels
title('Upper and Lower flows vs Pressure in the thorax')
xlabel('Pressure (dynes/cm^2)')
ylabel('Flow (cm^3/s)')
hold on
% line labels
label(p1,'case 1: Psv = 0','location','center',...
    'verticalalignment','bottom','fontweight','bold')
label(p2,'case 2:Psv = Pra','location','center',...
    'verticalalignment','bottom','fontweight','bold')
hold off

