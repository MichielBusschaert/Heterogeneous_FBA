function [ndf_idx,vxndf_idx] = getVariableIdxLHFBA(model)
    %Variable sizes
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    
    %Index ranges
    ndf_idx = 1:Ncells;
    vxndf_idx = reshape(Ncells+1:Ncells*(Nrxns+1),Nrxns,Ncells);
end