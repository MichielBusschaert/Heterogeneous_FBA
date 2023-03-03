function [Aeq,beq] = getConstraintBiomassProductionLHFBA(model,mu,B0)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    deltax = model.ub-model.lb;
    Aeq = zeros(1,Ncells*(Nrxns+1));
    beq = mu*B0;
    [~,vxndf_idx] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    dimS = length(size(model.S));
    Sp = [];
    for cell_idx = 1:Ncells
        %Assign stoichiometric matrix
        if dimS > 2
            Sp = model.S(model.sizeYmet+model.sizeXmet+1:model.sizeYmet+model.sizeXmet+model.sizePmet,:,cell_idx);
        elseif isempty(Sp)
            Sp = model.S(model.sizeYmet+model.sizeXmet+1:model.sizeYmet+model.sizeXmet+model.sizePmet,:);
        end
        Aeq(:,vxndf_idx(:,cell_idx)) = Sp*deltax(cell_idx);
    end
end