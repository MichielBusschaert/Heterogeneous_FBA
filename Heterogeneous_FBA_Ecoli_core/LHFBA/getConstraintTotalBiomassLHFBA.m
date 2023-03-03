function [Aeq,beq] = getConstraintTotalBiomassLHFBA(model,B0)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    Aeq = zeros(1,Ncells*(Nrxns+1));
    beq = B0;
    deltax = model.ub-model.lb;
    [ndf_idx,~] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    for cell_idx = 1:Ncells
        Aeq(:,ndf_idx(cell_idx)) = model.xcells(cell_idx)*deltax(cell_idx);
    end
end