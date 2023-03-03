function [Aeq,beq] = getConstraintMetabolismLHFBA(model)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    Aeq = zeros(Ncells*model.sizeXmet,Ncells*(Nrxns+1));
    beq = zeros(Ncells*model.sizeXmet,1);
    [~,vxndf_idx] = getVariableIdxLHFBA(model);
    
    %Iterate over cells
    dimS = length(size(model.S));
    Sx = [];        
    for cell_idx = 1:Ncells
        %Assign stoichiometric matrix
        if dimS > 2
            Sx = model.S(model.sizeYmet+1:model.sizeYmet+model.sizeXmet,:,cell_idx);
        elseif isempty(Sx)
            Sx = model.S(model.sizeYmet+1:model.sizeYmet+model.sizeXmet,:);
        end
        
        %Assign constraint
        cell_range = (cell_idx-1)*model.sizeXmet+1:cell_idx*model.sizeXmet;
        Aeq(cell_range,vxndf_idx(:,cell_idx)) = Sx;
    end
end