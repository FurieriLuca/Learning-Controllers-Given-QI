%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Learning distributed controllers given QI.
%% Luca Furieri, Yang Zheng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;



OutputCost_create_system_distributed;  %%sets up the system
OutputCost_convexprogram_outputfeedback;  %% solves the program with convex programming, so that we know the optimal cost and solution
optimal_cost = value(cost); %%%optimal cost




%ORACLE SETUP
OutputCost_create_cost;    %%sets up the symbolic cost. Not used later
eval_cost = matlabFunction(cost,'Vars',{vec_K}); %% This will be used to track evolution of the algorithm.



parameters = parameters_optimal - 1*ones(cardinality,1); %%initial parameters
parameters_initial = parameters;

fprintf('Starting the learning...from a cost of %d\n',eval_cost(parameters_initial))
pause(3)






rounds_max = 10;  %%we take the average of steps needed to achieve accuracy eps over rounds_max runs
epsilons = logspace(log10(0.2),log10(0.02),7) %% vector of logarithimically separated precisions



countnans = 0;  %%counter for divergent runs
countoutside = 0; %% counter for runs where the iterates exit G0

Taverages = zeros(1,size(epsilons,2)); %% This will contain the average number of iterates to achieve any precision over rounds_max runs of the algorithm

%eta_dividers=0.2*[1;1;1;1;1;1;1];   %%empirically found for good convergence at every considered precision
%r_dividers=1*[1;1;1;1;1;1;1];

eta = 0.0005; %%
r = sqrt(0.01);




running_window = 10000;   %% to monitor the progress in long runs. Will be used to pring max, average and min over a window of 3000 iterates
last_iterates_running = 1000000000*ones(1,running_window); %start from high value

i=1;  %initiate iteration

%%% SET THIS TO 1 IF YOU WANT TO PLOT COST; MAKES CODE MUCH SLOWER
plotting_cost = 1;
%%%

if(plotting_cost == 1)
    costItr = [];  %% will store the costs over time
    costRounds = zeros(1000000,rounds_max);
end
for(rounds = 1:rounds_max)
    count_success = 1; %%counts reaching a new precision level. Will go up to size(epsilons).
    costItr = [];
    fprintf('\n\nStarting round %d of precision %d\n',rounds,eps)
    while(true)
        
        %MONITORING OF PROGRESS OF REAL COST FUNCTION
        expected_cost_evaluated = eval_cost(parameters);
        if(plotting_cost == 1)
            costItr = [costItr;expected_cost_evaluated];
        end
        last_iterates_running = circshift(last_iterates_running,-1); %updates running window
        last_iterates_running(end) = expected_cost_evaluated;
        
        
        %SAMPLING
        [U, cost_sample] = sampling(cardinality, positions, r, parameters,n,m,p,N,x0a,x0b,wa,wb,va,vb,C_b,P11,P12,M_b,R_b); %%we sample the noisy cost to get cost_sample and U
        
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
    
    %%
    if(plotting_cost==1)
        costRounds([1:size(costItr,1)],rounds) = costItr;
    end
    
    last_iterates_running = 1000000000*ones(1,running_window);       % Reset running window
    Taverages = Taverages + Tsuccess; %%compute average steps over different rounds
    parameters = parameters_initial; %% resets parameters
end
Taverages = Taverages/rounds_max %% Compute actual average after finishing all rounds




%%PLOTS %%%First we handle the data for the plotting
final_steps = zeros(rounds_max,1);
iterate=1;
found = 0;
if(plotting_cost == 1)
    for(rounds=1:rounds_max)
        while(found == 0)
            if(costRounds(iterate,rounds) == 0)
                found = 1;
                final_steps(rounds) = iterate;
            end
            iterate = iterate + 1;
        end
        found = 0;
        iterate = 1;
    end
end

%%find max
maxes = zeros(rounds_max,1);
for(i=1:max(final_steps))
    maxes(i) = max(costRounds(i,:));
end

%find min
for(i=1:min(final_steps))
    mins(i) = min(costRounds(i,:));
end
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

extended_y = [mins([1:min(final_steps)-1])';optimal_cost+0.019*ones(max(final_steps)-min(final_steps),1)];
extended_y = extended_y - 0.001*ones(size(extended_y,1),1);


%%On a second figure, we plot the iterates for one run.
if(plotting_cost == 1)
    figure(2)
    plot(maxes([1:max(final_steps)-1]),'-k','linewidth',1.5);
    hold on
    grid on
    plot(mins([1:min(final_steps)-1]),'-k','linewidth',1.5);
    plot(linspace(1,max(final_steps)),optimal_cost*ones(100,1),'--r','linewidth',2)
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(1))*ones(100,1),':g','linewidth',1.5)
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(2))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(3))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(4))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(5))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(6))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    plot(linspace(1,max(final_steps)),(optimal_cost+epsilons(7))*ones(100,1),':g','linewidth',1.5,'HandleVisibility','off')
    axis([0 max(final_steps) optimal_cost-0.001 eval_cost(parameters_initial)+0.001])
    
%     
   % x = [[1:max(final_steps)-1]';[1:max(final_steps)-1]'];
   % y = [extended_y;maxes([1:max(final_steps)-1])];
   % fill(x,y,'g')
    
end

save('DATA')
%grid on
%scatter(epsilons,Ts,50,'bo','Linewidth',2)
%set(gca, 'xscale', 'log', 'yscale', 'log');




%fprintf('DONE!\n')
%final_cost=eval_cost(parameters)
%initial_cost=eval_cost(parameters_initial)





