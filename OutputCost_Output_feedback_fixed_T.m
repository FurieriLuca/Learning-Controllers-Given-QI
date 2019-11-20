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





parameters = parameters_optimal + 1*ones(cardinality,1);
parameters_initial = parameters; %%% chooses initial parameters%[1.3828;-1.6278;1.4522]+0.5*[1;1;1];%1*[5;-5;0]%[1.38+1;-1.6;1.40];%10*(rand(cardinality,1)-rand(cardinality,1));
fprintf('Starting the learning...from a cost of %f\n',eval_cost(parameters_initial))
pause(3)






rounds_max = 2;  %%we take the average of steps needed to achieve accuracy eps over rounds_max runs

epsilons = round(logspace(log10(0.5),log10(10),10),1) %% creates a logarithmic space of precisions
for(i = 1:size(epsilons,2))
    epsilons(i) = 1/epsilons(i);
end


parameters_success = zeros(cardinality,rounds_max,size(epsilons,2));



countnans = 0;  %%counter for divergent runs
countoutside = 0; %% counter for runs where the iterates exit G0
count_succeeded = zeros(size(epsilons,2),1);
avg_final_cost = zeros(size(epsilons,2),1);


epsilons = logspace(log10(0.2),log10(0.01),7) %% space of precisions
eta_dividers=[0.02;0.02;0.02;0.02;0.02;0.02;0.02];   %%empirically found for good convergence at every considered precision
r_dividers=[0.1;0.1;0.1;0.1;0.1;0.1;0.1]; 
for(precisions = 1:size(epsilons,2))
    
    eps = epsilons(precisions);  %chooses the precision
    fprintf('\n\n\n***NEW PRECISION*** : %f\n',eps)
    
    eta = eps^2/eta_dividers(precisions); %115% %% stepsize (goes with eps^2)
    r = sqrt(eps)/r_dividers(precisions);   %% smoothing radius (goes with sqrt(eps) )
    samples_number = 1;  %%mini-batch size
    
    %% eta2=4*0.000115, eta1.5= . , eta1=0.000115, eta0.5=0.000115/2
    
    
    T=floor((4/eta*log(120*(eval_cost(parameters_initial)-optimal_cost)/eps))*367);
    
    window = 2000;  %% to monitor the progress in long runs, we will print the average of the real cost over the past "window" iterates
    running_window = window;
    
    last_iterates_running = 1000000000*ones(1,running_window);
    
    i=1;
    for(rounds = 1:rounds_max)
        fprintf('\n\nStarting round %f of precision %f, running for %d steps\n',rounds,eps,T)
        while(i<=T)
            
            %MONITORING OF PROGRESS OF REAL COST FUNCTION
            expected_cost_evaluated = eval_cost(parameters);
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
                fprintf('RE-Starting round %f of precision %f\n',rounds,eps)
                i = 1;
                parameters = parameters_initial;
                last_iterates_running = 1000000000*ones(1,running_window);
            elseif((eval_cost(parameters)-optimal_cost)>10*(eval_cost(parameters_initial)-optimal_cost)) %%if it goes outside G0... restart the round
                countoutside = countoutside+1;
                fprintf('!!!!!!Went outside G0!!!!!!\n')
                fprintf('RE-Starting round %f of precision %f\n',rounds,eps)
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
            
            i = i+1;
        end
        
        if(eval_cost(parameters)-optimal_cost<eps)
            final_cost = eval_cost(parameters);
            avg_final_cost(precisions) = avg_final_cost(precisions) + final_cost;
            fprintf('!Succeeded! final cost %f\n',eval_cost(parameters))
            count_succeeded(precisions)=count_succeeded(precisions)+1;

        else
            fprintf('!FAILED! final cost %f\n',eval_cost(parameters))
        end
        i = 1;
        last_iterates_running = 1000000000*ones(1,running_window);        %variance = 0;
        parameters = parameters_initial;
    end
    avg_final_cost(precisions) = avg_final_cost(precisions)/rounds_max
    fprintf('%d/%d attempts were succesful for precision eps = %f . \n',count_succeeded(precisions),rounds_max,epsilons(precisions));
    
end


Ts = zeros(1,size(epsilons,2));
for(i=1:size(epsilons,2))
    eta = epsilons(i)^2*0.0125
    Ts(i) = floor((4/eta*log(120*(eval_cost(parameters_initial)-optimal_cost)/epsilons(i)))*1.8080)
end


figure(1)
plot(epsilons,Ts,'o-','linewidth',2,'linestyle','--','markersize',10,'markerfacecolor','c')
grid on
%scatter(epsilons,Ts,50,'bo','Linewidth',2)
set(gca, 'xscale', 'log', 'yscale', 'log');
xlim([0.01 1]) ;figure(1)
ylim([0.5*Ts(1) 2*Ts(end)]) ;
%Coeffs = polyfit(log(epsilons([1:end])),log(Ts([1:end])),2);
%hold on
%grid on
%loglog(epsilons, exp(polyval(Coeffs, log(epsilons))),'LineStyle','--','LineWidth',1);





%fprintf('DONE!\n')
%final_cost=eval_cost(parameters)
%initial_cost=eval_cost(parameters_initial)





