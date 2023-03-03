function [Aeq,beq] = getConstraintHeterogeneityLHFBA(model,mu)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    Aeq = zeros(Ncells,Ncells*(Nrxns+1));
    beq = zeros(Ncells,1);
    deltax = model.ub-model.lb;
    division_gamma_ct = model.division_gamma_ct;
    division_gamma_mu = model.division_gamma_mu;
    division_beta_mat = model.division_beta_mat;
    [ndf_idx,vxndf_idx] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    for cell_idx = 1:Ncells
        xcell = model.xcells(cell_idx);
        
        %Balanced growth term (quasi-steady-state)
        Aeq(cell_idx,ndf_idx(cell_idx)) = mu;
        
        %Growth term - 1st order upwind scheme: dfdx ~~ (f{i}-f{i-1})/dx 
        Aeq(cell_idx,vxndf_idx(:,cell_idx)) = model.c'/deltax(cell_idx);
        if cell_idx ~= 1
            Aeq(cell_idx,vxndf_idx(:,cell_idx-1)) = -model.c'/deltax(cell_idx);
        end
        
%         %Growth term - 2nd order upwind scheme: dfdx ~~
%         %(3*f{i}-4*f{i-1}+f{i-2})/(2*dx)
%         Aeq(cell_idx,vxndf_idx(:,cell_idx)) = 1.5*model.c'/deltax(cell_idx);
%         if cell_idx > 1
%             Aeq(cell_idx,vxndf_idx(:,cell_idx-1)) = -2*model.c'/deltax(cell_idx);
%         end
%         if cell_idx > 2
%             Aeq(cell_idx,vxndf_idx(:,cell_idx-2)) = .5*model.c'/deltax(cell_idx);
%         end

%         %Growth term - hybrid upwind
%         if cell_idx > 2
%             Aeq(cell_idx,vxndf_idx(:,cell_idx)) = 1.5*model.c'/deltax(cell_idx);
%             Aeq(cell_idx,vxndf_idx(:,cell_idx-1)) = -2*model.c'/deltax(cell_idx);
%             Aeq(cell_idx,vxndf_idx(:,cell_idx-2)) = .5*model.c'/deltax(cell_idx);
%         elseif cell_idx > 1
%             Aeq(cell_idx,vxndf_idx(:,cell_idx)) = model.c'/deltax(cell_idx);
%             Aeq(cell_idx,vxndf_idx(:,cell_idx-1)) = -model.c'/deltax(cell_idx);
%         else
%             Aeq(cell_idx,vxndf_idx(:,cell_idx)) = model.c'/deltax(cell_idx);
%         end
        
        %Division term - death
        Aeq(cell_idx,ndf_idx(cell_idx)) = Aeq(cell_idx,ndf_idx(cell_idx)) + division_gamma_ct(xcell);
        Aeq(cell_idx,vxndf_idx(:,cell_idx)) = Aeq(cell_idx,vxndf_idx(:,cell_idx)) + model.c'*division_gamma_mu(xcell);
        
        %Division term - birth
        Aeq(cell_idx,ndf_idx(cell_idx)) = Aeq(cell_idx,ndf_idx(cell_idx)) - division_gamma_ct(xcell)*division_beta_mat(cell_idx,cell_idx)*deltax(cell_idx);
        Aeq(cell_idx,vxndf_idx(:,cell_idx)) = Aeq(cell_idx,vxndf_idx(:,cell_idx)) - model.c'*division_gamma_mu(xcell)*division_beta_mat(cell_idx,cell_idx)*deltax(cell_idx);
        for break_idx = cell_idx+1:Ncells
            xbreak = model.xcells(break_idx);
            Aeq(cell_idx,ndf_idx(break_idx)) = -division_gamma_ct(xbreak)*division_beta_mat(cell_idx,break_idx)*deltax(break_idx);
            Aeq(cell_idx,vxndf_idx(:,break_idx)) = -model.c'*division_gamma_mu(xbreak)*division_beta_mat(cell_idx,break_idx)*deltax(break_idx);
        end
    end
end