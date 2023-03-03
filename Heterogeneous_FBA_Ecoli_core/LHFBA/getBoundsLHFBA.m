function [lb,ub] = getBoundsLHFBA(model)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    lb = -Inf*ones(Ncells*(Nrxns+1),1);
    ub = Inf*ones(Ncells*(Nrxns+1),1);
    [ndf_idx,vxndf_idx] = getVariableIdxLHFBA(model);

    %NDF positivity
    lb(ndf_idx) = zeros(Ncells,1);
    
    %Reaction irreversibility
    for k=1:Nrxns
        if (model.rev(k)==0)
            lb(vxndf_idx(k,:)) = 0;
        end
    end
end