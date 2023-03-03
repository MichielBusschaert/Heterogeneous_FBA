function [Aineq,bineq] = getConstraintExchangeLHFBA(model)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    deltax = model.ub-model.lb;
    num_idx = ~(isinf(model.uptakeLim) | isnan(model.uptakeLim));
    Aineq = zeros(sum(num_idx),Ncells*(Nrxns+1));
    bineq = -model.uptakeLim(num_idx);
    [~,vxndf_idx] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    dimS = length(size(model.S));
    Sy = [];
    for cell_idx = 1:Ncells
        %Assign stoichiometric matrix
        if dimS > 2
            Sy = model.S(num_idx,:,cell_idx);
        elseif isempty(Sy)
            Sy = model.S(num_idx,:);
        end
        Aineq(:,vxndf_idx(:,cell_idx)) = -Sy*deltax(cell_idx);
    end
end