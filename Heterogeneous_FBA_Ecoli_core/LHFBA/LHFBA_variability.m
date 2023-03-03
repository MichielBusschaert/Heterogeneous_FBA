%% LHFBA_variability
function [ndf_lb,v_lb,ndf_ub,v_ub] = LHFBA_variability(model,mu,B0,hfba_options,lpsolver)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + 1;
    [ndf_lb,ndf_ub] = deal(zeros(Ncells,1));
    [v_lb,v_ub] = deal(zeros(Nrxns,Ncells));
    [ndf_idx,vxndf_idx] = getVariableIdxLHFBA(model);

    %Formulate constraints
    [Aineq,bineq,Aeq,beq,lb,ub] = formulateConstraintsLHFBA(model,mu,B0,hfba_options);
    
    %Formulate objective function
    f = formulateObjectiveLHFBA(model);
    
    %Solver options
    switch lpsolver
        case 'matlab'
            disp('Solving using linprog');
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
            elseif lpflag < 0 % Infeasible problem
                disp('Infeasible LHFBA problem, consider using other initial points.')
                ndf = [];
                vxndf = [];
            else
                disp('LHFBA not solved within maximum number of iterations')
                ndf = [];
                vxndf = [];
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

            % create cplex object
            lpProb_fl = Cplex;
            % set parameters
            lpProb_fl.Param.simplex.tolerances.optimality.Cur = 1e-18;
            lpProb_fl.Param.simplex.tolerances.feasibility.Cur = 1e-18;
            lpProb_fl.Param.feasopt.tolerance.Cur = 1e-18;
            lpProb_fl.Param.sifting.display.Cur = 0;
            lpProb_fl.Param.emphasis.numerical.Cur = 1;
            lpProb_fl.Param.simplex.display.Cur = 0;
            lpProb_fl.Param.tune.display.Cur = 3;
            lpProb_fl.Param.barrier.convergetol.Cur = 1e-20;
            lpProb_fl.Param.barrier.display.Cur = 0;
            lpProb_fl.Param.threads.Cur = 12;
            lpProb_fl.Param.workmem.Cur = 2048;     
            lpProb_fl.Param.read.scale.Cur = 0;
            lpProb_fl.DisplayFunc =[];

            % add columns and rows: NDF
            lpProb.addCols(f,[],lb,ub);
            lpProb.addRows([beq;-Inf(size(bineq))],[Aeq;Aineq],[beq;bineq]);

            % add columns and rows: NDF
            Aineq_fl = [[Aineq,-bineq];
                [-eye(Ncells*(Nrxns+1)),lb];
                [eye(Ncells*(Nrxns+1)),-ub]];
            Aineq_fl(any(isinf(Aineq_fl),2),:) = [];
            Aeq_fl = [Aeq,-beq];
            bineq_fl = zeros(size(Aineq_fl,1),1);
            beq_fl = zeros(size(Aeq_fl,1),1);
            lb_fl = [-Inf(Ncells*(Nrxns+1),1);0];
            ub_fl = Inf(Ncells*(Nrxns+1)+1,1);
            f_fl = [f;0];
            lpProb_fl.addCols(f_fl,[],lb_fl,ub_fl);
            lpProb_fl.addRows([beq_fl;-Inf(size(bineq_fl))],[Aeq_fl;Aineq_fl],[beq_fl;bineq_fl]);

            %Iterate over cells
            for cell_idx = 1:Ncells
                %NDF variability
                f = zeros(Ncells*(Nrxns+1),1);
                f(ndf_idx(cell_idx)) = 1;
                lpProb.Model.obj = f;
                
                %Solve min
                lpProb.Model.sense = 'minimize';
                lpProb.solve();
                lp_solution = lpProb.Solution.x;
                ndf = lp_solution(1:Ncells);
                ndf_lb(cell_idx) = ndf(cell_idx);
                
                %Solve max
                lpProb.Model.sense = 'maximize';
                lpProb.solve();
                lp_solution = lpProb.Solution.x;
                ndf = lp_solution(1:Ncells);
                ndf_ub(cell_idx) = ndf(cell_idx);

                %Flux variability
                d_fl = zeros(1,Ncells*(Nrxns+1)+1);
                d_fl(ndf_idx(cell_idx)) = model.xcells(cell_idx);
                lpProb_fl.addRows(1,d_fl,1);
                for flux_idx = 1:Nrxns
                    f_fl = zeros(Ncells*(Nrxns+1)+1,1);
                    f_fl(vxndf_idx(flux_idx,cell_idx)) = 1;
                    lpProb_fl.Model.obj = f_fl;

                    %Solve min
                    lpProb_fl.Model.sense = 'minimize';
                    lpProb_fl.solve();
                    lp_solution = lpProb_fl.Solution.x;
                    rndf = lp_solution(1:Ncells);
                    rvxn = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
                    v_lb(flux_idx,cell_idx) = rvxn(flux_idx,cell_idx)/rndf(cell_idx)/model.xcells(cell_idx);

                    %Solve max
                    lpProb_fl.Model.sense = 'maximize';
                    lpProb_fl.solve();
                    lp_solution = lpProb_fl.Solution.x;
                    rndf = lp_solution(1:Ncells);
                    rvxn = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
                    v_ub(flux_idx,cell_idx) = rvxn(flux_idx,cell_idx)/rndf(cell_idx)/model.xcells(cell_idx);

%                     f_fl = zeros(Ncells*(Nrxns+1)+1,1);
%                     f_fl(vxndf_idx(flux_idx,cell_idx)) = 1;
%                     lpProb_fl.Model.obj = f_fl;
% 
%                     %Solve min
%                     lpProb.Model.sense = 'minimize';
%                     lpProb.solve();
%                     lp_solution = lpProb.Solution.x;
%                     ndf = lp_solution(1:Ncells);
%                     vxndf = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
%                     v = vxndf./repmat(ndf.*model.xcells,1,model.sizeYrxn+model.sizeXrxn+1)';
%                     v_lb(flux_idx,cell_idx) = v(flux_idx,cell_idx);
% 
%                     %Solve max
%                     lpProb.Model.sense = 'maximize';
%                     lpProb.solve();
%                     lp_solution = lpProb.Solution.x;
%                     ndf = lp_solution(1:Ncells);
%                     vxndf = reshape(lp_solution(Ncells+1:Ncells*(Nrxns+1)),Nrxns,Ncells);
%                     v = vxndf./repmat(ndf.*model.xcells,1,model.sizeYrxn+model.sizeXrxn+1)';
%                     v_ub(flux_idx,cell_idx) = v(flux_idx,cell_idx);
                end
                lpProb_fl.delRows(size(lpProb_fl.Model.A,1));
            end
    end
end 