%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Learning distributed controllers given QI.
%% Luca Furieri, Yang Zheng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;

OutputCost_convexprogram_outputfeedback;  %% solves the program with convex programming, so that we know the optimal cost and solution
optimal_cost = value(cost); %%%optimal cost

%{
figure(2)
hold on
scatter3(K_opt(1,1),K_opt(2,1),K_opt(2,4),'kx')
xlim([-5 5]);
ylim([-5 5]);
zlim([-5 5]);
%}


%ORACLE SETUP
OutputCost_create_cost;    %%sets up the symbolic cost. Not used later
eval_cost = matlabFunction(cost,'Vars',{vec_K}); %%



parameters = parameters_optimal - 1*ones(cardinality,1);
parameters_initial = parameters; %%% chooses initial parameters%[1.3828;-1.6278;1.4522]+0.5*[1;1;1];%1*[5;-5;0]%[1.38+1;-1.6;1.40];%10*(rand(cardinality,1)-rand(cardinality,1));
fprintf('Starting the learning...from a cost of %d\n',eval_cost(parameters_initial))
pause(3)






rounds_max = 1;  %%we take the average of steps needed to achieve accuracy eps over rounds_max runs
epsilons = logspace(log10(0.2),log10(0.02),7) %% space of precisions



parameters_success = zeros(cardinality,rounds_max,size(epsilons,2));



countnans = 0;  %%counter for divergent runs
countoutside = 0; %% counter for runs where the iterates exit G0

Taverages = zeros(1,size(epsilons,2));

eta_dividers=0.2*[1;1;1;1;1;1;1];   %%empirically found for good convergence at every considered precision
r_dividers=1*[1;1;1;1;1;1;1];
%for(precisions = 1:size(epsilons,2))

%eps = epsilons(precisions);  %chooses the precision
%fprintf('\n\n\n***NEW PRECISION*** : %d\n',eps)

eta = 0.01^2/eta_dividers(1); %115% %% stepsize (goes with eps^2)
r = sqrt(0.01)/r_dividers(1);   %% smoothing radius (goes with sqrt(eps) )
samples_number = 1;  %%mini-batch size

%% eta2=4*0.000115, eta1.5= . , eta1=0.000115, eta0.5=0.000115/2



window = 3000;  %% to monitor the progress in long runs, we will print the average of the real cost over the past "window" iterates
running_window = window;

last_iterates_running = 1000000000*ones(1,running_window);

%Taverage = 0;
i=1;
costItr = [];
for(rounds = 1:rounds_max)
    count_success=1;
    fprintf('\n\nStarting round %d of precision %d\n',rounds,eps)
    while(true)
        
        %MONITORING OF PROGRESS OF REAL COST FUNCTION
        expected_cost_evaluated = eval_cost(parameters);
        costItr = [costItr;expected_cost_evaluated];
        last_iterates_running = circshift(last_iterates_running,-1);
        last_iterates_running(end) = expected_cost_evaluated;
        
        
        %SAMPLING
        grad_estimate = zeros(cardinality,1);
        for(samples = 1:samples_number)
            OutputCost_sample_cost; %%we sample the noisy cost to get cost_sample and U
            grad_estimate = grad_estimate+cost_sample*U;  %gradient estimate
        end
        grad_estimate = grad_estimate*cardinality/r^2/samples_number; %correctly scaled gradient estimate
        
        
        %%STEP
        parameters = parameters-eta*grad_estimate;
        
        
        %%CHECKS DIVERGENCE AND GOING OUT OF G0
        if(isnan(parameters)~=zeros(cardinality,1))  %%if diverges... restart the round
            countnans = countnans+1;
            fprintf('!!!!NAN!!! reset the round\n')
            fprintf('RE-Starting round %d of precision %d\n',rounds,epsilons(count_success))
            i = 1;
            parameters = parameters_initial;
            last_iterates_running = 1000000000*ones(1,running_window);
        elseif((eval_cost(parameters)-optimal_cost)>10*(eval_cost(parameters_initial)-optimal_cost)) %%if it goes outside G0... restart the round
            countoutside = countoutside+1;
            fprintf('!!!!!!Went outside G0!!!!!!\n')
            fprintf('RE-Starting round %d of precision %d\n',rounds,epsilons(count_success))
            i = 1;
            parameters = parameters_initial;
            last_iterates_running = 1000000000*ones(1,running_window);
        end
        %}
        
        
        
        %%PRINTING
        if(mod(i,window)==0) %%%print the progress every window-th step
            max_mean_min_window = [max(last_iterates_running), mean(last_iterates_running), min(last_iterates_running)]
            i
            
            %real_expected_cost = 0;
        end
        
        %STOPPING CRITERION
        if(expected_cost_evaluated-optimal_cost<epsilons(count_success))
            
            
            fprintf('\n\n\n***entered precision bound %d at step %d ***\n\n\n',epsilons(count_success),i)
            Tsuccess(count_success) = i;
            count_success = count_success + 1;
            if(count_success == size(epsilons,2)+1)
                break;
            end
            %parameters_success(:,rounds,precisions) = parameters;
            %plot_parameter
            %break;
        end
        
        i = i+1;
    end
    i = 1;
    
    %avg = 0;
    last_iterates_running = 1000000000*ones(1,running_window);        %variance = 0;
    Taverages = Taverages + Tsuccess;
    parameters = parameters_initial;
end
Taverages = Taverages/rounds_max %% average steps needed over rounds_max runs
%Taverage = 0;
%end

%openfig('figures/T_fixed_AND_firststop')
figure(1)
hold on
grid on
%plot(epsilons,Ts,'ro--','linewidth',2,'markersize',10,'markerfacecolor','m')

plot(epsilons,Taverages,'bo-','linewidth',2,'linestyle','--','markersize',15,'markerfacecolor','b')
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([0.02 0.2]);
ylim([65000 11057335]) ;



%Coeffs = polyfit(log(epsilons),log(Taverages),1);

%loglog(epsilons, exp(polyval(Coeffs, log(epsilons))),'LineStyle','-','LineWidth',1);


Ts = zeros(1,size(epsilons,2));
for(i=1:size(epsilons,2))
    eta = epsilons(i)^2*0.0125
    Ts(i) = floor((4/eta*log(120*(eval_cost(parameters_initial)-optimal_cost)/epsilons(i)))*1.8080)
end
plot(epsilons,Ts,'ro-','linewidth',2,'linestyle','--','markersize',15,'markerfacecolor','r')

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



%grid on
%scatter(epsilons,Ts,50,'bo','Linewidth',2)
%set(gca, 'xscale', 'log', 'yscale', 'log');




%fprintf('DONE!\n')
%final_cost=eval_cost(parameters)
%initial_cost=eval_cost(parameters_initial)





