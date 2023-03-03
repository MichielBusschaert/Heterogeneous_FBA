function [ndf,vxndf,sol_object] = LHFBA_function(model,mu,B0,hfba_options,lpsolver)
    %Formulate constraints
    [Aineq,bineq,Aeq,beq,lb,ub] = formulateConstraintsLHFBA(model,mu,B0,hfba_options);
    
    %Formulate objective function
    f = formulateObjectiveLHFBA(model);
    
    %Solver options
    switch lpsolver
        case 'matlab'
            disp('Solving using linprog');
            %options = optimoptions('fmincon','ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',6e5,'MaxIterations',2e5,'BarrierParamUpdate','predictor-corrector');

            %Solve non-linear program

            %[lp_solution,~,lpflag] = fmincon(f,x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
            options = optimoptions('linprog','OptimalityTolerance',1e-10,'ConstraintTolerance',1e-6);
            [lp_solution,~,lpflag] = linprog(-f,Aineq,bineq,Aeq,beq,lb,ub,options);
            if lpflag > 0 % Solution found
                if isempty(find(lp_solution)) % Check if it is not the zero flux solution
                    disp('Trivial LHFBA solution obtained')
                else
                    disp('Non-trivial LHFBA solution obtained')
                end
                Ncells = model.sizeCells;
                Nrxns = model.sizeYrxn + model.sizeXrxn + 1;

                ndf = lp_solution(1:Ncells);
                vxndf = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
                sol_object = lpflag;
            elseif lpflag < 0 % Infeasible problem
                disp('Infeasible LHFBA problem, consider using other initial points.')
                ndf = [];
                vxndf = [];
                sol_object = lpflag;
            else
                disp('LHFBA not solved within maximum number of iterations')
                ndf = [];
                vxndf = [];
                sol_object = lpflag;
            end
        case 'cplex'
            disp('Solving using cplex');
            % create cplex object
            lpProb = Cplex;
            % set parameters
            lpProb.Param.simplex.tolerances.optimality.Cur = 1e-18;
            lpProb.Param.simplex.tolerances.feasibility.Cur = 1e-18;
            lpProb.Param.feasopt.tolerance.Cur = 1e-18;
            lpProb.Param.sifting.display.Cur = 0;
            lpProb.Param.emphasis.numerical.Cur = 1;
            lpProb.Param.simplex.display.Cur = 0;
            lpProb.Param.tune.display.Cur = 3;
            lpProb.Param.barrier.convergetol.Cur = 1e-20;
            lpProb.Param.barrier.display.Cur = 0;
            lpProb.Param.threads.Cur = 12;
            lpProb.Param.workmem.Cur = 2048;     
            lpProb.Param.read.scale.Cur = 0;
            lpProb.DisplayFunc =[];
            
            % add bounds and objective
            lpProb.addCols(f,[],lb,ub);
            lpProb.Model.sense = 'maximize';
            lpProb.addRows([beq;-Inf(size(bineq))],[Aeq;Aineq],[beq;bineq]);
            
            % solve LP
            lpProb.solve();
            if lpProb.Solution.status == 1 % Solution found
                lp_solution = lpProb.Solution.x;
                if isempty(find(lp_solution)) % Check if it is not the zero flux solution
                    disp('Trivial LHFBA solution obtained')
                else
                    disp('Non-trivial LHFBA solution obtained')
                end
                Ncells = model.sizeCells;
                Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;

                ndf = lp_solution(1:Ncells);
                vxndf = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
                sol_object = lpProb;
            elseif lpProb.Solution.status == 3 % Infeasible problem
                disp('Infeasible HFBA problem, consider using other initial points.')
                ndf = [];
                vxndf = [];
                sol_object = lpProb;
            end
    end
end