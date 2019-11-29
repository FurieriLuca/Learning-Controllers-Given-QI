%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Learning distributed controllers given QI.
%% Luca Furieri, Yang Zheng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;



OutputCost_create_system_distributed_ACC;  %%sets up the system





%ORACLE SETUP
OutputCost_create_cost;    %%sets up the symbolic cost. Not used later
%%%optimal cost
parameters_optimal = [0.275208;1.13539];

eval_cost = matlabFunction(cost,'Vars',{[a;b]}); %% This will be used to track evolution of the algorithm.
optimal_cost=eval_cost(parameters_optimal);
cardinality=2;


parameters = parameters_optimal + 0.5*ones(cardinality,1); %%initial parameters
parameters_initial = parameters;

fprintf('Starting the learning...from a cost of %d\n',eval_cost(parameters_initial))
pause(3)






rounds_max = 5;  %%we take the average of steps needed to achieve accuracy eps over rounds_max runs
epsilons = logspace(log10(0.002),log10(0.0002),7) %% vector of logarithimically separated precisions



countnans = 0;  %%counter for divergent runs
countoutside = 0; %% counter for runs where the iterates exit G0

Taverages = zeros(1,size(epsilons,2)); %% This will contain the average number of iterates to achieve any precision over rounds_max runs of the algorithm

%eta_dividers=0.2*[1;1;1;1;1;1;1];   %%empirically found for good convergence at every considered precision
%r_dividers=1*[1;1;1;1;1;1;1];

eta = 0.00000005; %%
r = sqrt(0.01)*0.89;




running_window = 10000;   %% to monitor the progress in long runs. Will be used to pring max, average and min over a window of 3000 iterates
last_iterates_running = 1000000000*ones(1,running_window); %start from high value

i=1;  %initiate iteration

%%% SET THIS TO 1 IF YOU WANT TO PLOT COST; MAKES CODE MUCH SLOWER
plotting_cost = 0; 
%%%

if(plotting_cost == 1)
    costItr = [];  %% will store the costs over time
end
for(rounds = 1:rounds_max)
    count_success = 1; %%counts reaching a new precision level. Will go up to size(epsilons).
    costItr = [];
    fprintf('\n\nStarting round %d of precision %d\n',rounds,epsilons(count_success))
    while(true)
        
        %MONITORING OF PROGRESS OF REAL COST FUNCTION
        expected_cost_evaluated = eval_cost(parameters);
        if(plotting_cost == 1)
            costItr = [costItr;expected_cost_evaluated];
        end
        last_iterates_running = circshift(last_iterates_running,-1); %updates running window
        last_iterates_running(end) = expected_cost_evaluated;
        
        
        %SAMPLING
        [U, cost_sample] = sampling_ACC(cardinality, r, parameters,n,m,p,N,x0a,x0b,wa,wb,va,vb,C_b,P11,P12,M_b,R_b); %%we sample the noisy cost to get cost_sample and U
        
        %Gradient estimation
        grad_estimate = zeros(cardinality,1);
        grad_estimate = grad_estimate+cost_sample*U;  %gradient estimate
        grad_estimate = grad_estimate*cardinality/r^2; %correctly scaled gradient estimate
        
        
        %%STEP
        parameters = parameters-eta*grad_estimate;
        
        
        %%CHECKS GOING OUT OF G0
        if((eval_cost(parameters)-optimal_cost)>10*(eval_cost(parameters_initial)-optimal_cost)) %%if it goes outside G0... can consider a value different from "10"
            countoutside = countoutside+1;
            fprintf('!!!!!!Went outside G0!!!!!!\n')
            fprintf('RE-Starting round %d of precision %d\n',rounds,epsilons(count_success))
            i = 1;
            parameters = parameters_initial;
            last_iterates_running = 1000000000*ones(1,running_window);
            if(plotting_cost == 1)
                costItr = [];
            end
        end
        
        
        
        
        %%PRINTING
        if(mod(i,running_window)==0) %%%prints the progress every running_window-th step
            max_mean_min_window = [max(last_iterates_running), mean(last_iterates_running), min(last_iterates_running)]
            i
        end
        
        %STOPPING CRITERION
        if(expected_cost_evaluated-optimal_cost<epsilons(count_success))
            fprintf('\n\n\n***entered precision bound %d at step %d ***\n\n\n',epsilons(count_success),i)
            Tsuccess(count_success) = i;
            count_success = count_success + 1;  %we have reached a new precision level, continue iterations to reach the next precision level
            if(count_success == size(epsilons,2)+1)
                break;  %We have reached all precision levels. We can start a new round.
            end
        end
        
        i = i+1;
    end
    i = 1;
    
    last_iterates_running = 1000000000*ones(1,running_window);       % Reset running window
    Taverages = Taverages + Tsuccess; %%compute average steps over different rounds
    parameters = parameters_initial; %% resets parameters
end
Taverages = Taverages/rounds_max %% Compute actual average after finishing all rounds


%%PLOTS
%%First, we plot in log-log scale how many steps it took to reach each
%%precision level

figure(1)
hold on
grid on
plot(epsilons,Taverages,'bo-','linewidth',2,'linestyle','--','markersize',15,'markerfacecolor','b')
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([0.02 0.2]);
ylim([65000 11057335]) ;

% On the same picture, we plot the scaling with epsilon predicted by the
% theorem

Ts = zeros(1,size(epsilons,2));
for(i=1:size(epsilons,2))
    eta = epsilons(i)^2*0.0125
    Ts(i) = floor((4/eta*log(120*(eval_cost(parameters_initial)-optimal_cost)/epsilons(i)))*1.8080) %% The constants are chosen to make the first precision level correspond.
end
plot(epsilons,Ts,'ro-','linewidth',2,'linestyle','--','markersize',15,'markerfacecolor','r')


%%On a second figure, we plot the iterates for one run.
if(plotting_cost == 1)
    
    figure(2)
    plot(costItr,'-k','linewidth',1.5);
    hold on
    grid on
    plot(linspace(1,size(costItr,1)),optimal_cost*ones(100,1),'--r','linewidth',2)
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(1))*ones(100,1),':g','linewidth',1.5)
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(2))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(3))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(4))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(5))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(6))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,size(costItr,1)),(optimal_cost+epsilons(7))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
end


%grid on
%scatter(epsilons,Ts,50,'bo','Linewidth',2)
%set(gca, 'xscale', 'log', 'yscale', 'log');




%fprintf('DONE!\n')
%final_cost=eval_cost(parameters)
%initial_cost=eval_cost(parameters_initial)





