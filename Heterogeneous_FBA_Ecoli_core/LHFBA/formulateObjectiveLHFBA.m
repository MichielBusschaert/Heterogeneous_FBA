function f = formulateObjectiveLHFBA(model)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    f = zeros(Ncells*(Nrxns+1),1);
    deltax = model.ub-model.lb;
    [~,vxndf_idx] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    for cell_idx = 1:Ncells
        f(vxndf_idx(:,cell_idx)) = model.c'*deltax(cell_idx);
    end
end