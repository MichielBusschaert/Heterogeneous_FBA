function [ndf_sol,v_sol,mu_avg_sol,sol_object] = LHFBA_distribution(model,mu_avg_init,biomass,hfba_options)
    %Iterate with increasing growth rate until problem becomes infeasible
    solved_flag = false;
    ndf_sol = [];
    vxndf_sol = [];
    mu_avg_sol = 0;
    per_incr = 1;
    per_step = 10;
    per_min = 1e-4;

    while solved_flag == false
        disp(['Solving with mu = ',num2str(mu_avg_init),'.']);
        [ndf,vxndf,sol_object_tmp] = LHFBA_function(model,mu_avg_init,biomass,hfba_options,'cplex');
        
        if ~isempty(ndf)
            %Solution found but maybe not optimal growth rate
            ndf_sol = ndf;
            vxndf_sol = vxndf;
            mu_avg_sol = mu_avg_init;
            mu_avg_init = (1+per_incr)*mu_avg_init;
            sol_object = sol_object_tmp;
        else
            if ~isempty(ndf_sol)
                %Optimal growth rate exceded
                %Narrow growth rate if needed
                if per_incr >= per_min
                    mu_avg_init = mu_avg_init/(1+per_incr);
                    per_incr = per_incr/per_step;
                    mu_avg_init = (1+per_incr)*mu_avg_init;
%                 elseif any(ndf_sol<0)
%                     mu_avg_init = mu_avg_init/(1+per_incr);
                else
                    solved_flag = true;
                end
            else
                if mu_avg_init >= per_min
                    %Optimal growth rate exceded, but no prior solution
                    mu_avg_init = mu_avg_init/(1+per_incr);
                else
                    %No solution found --> mu_avg_init = 0
                    ndf_sol = [];
                    vxndf_sol = [];
                    mu_avg_sol = 0;
                    sol_object = sol_object_tmp;
                    solved_flag = true;
                end
            end
        end
    end
    
    %Recalculate optimal solution
%     [ndf_sol,vxndf_sol] = LHFBA_minimal_exchange(model,mu_avg_sol,biomass,hfba_options,'cplex');
    
    %Calculate fluxes
    if mu_avg_sol > 0
        v_sol = vxndf_sol./repmat(ndf_sol.*model.xcells,1,model.sizeYrxn+model.sizeXrxn+model.sizePrxn)';
    else
        v_sol = [];
    end
end