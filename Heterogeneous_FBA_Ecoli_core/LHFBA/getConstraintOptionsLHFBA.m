function [Aineq,bineq,Aeq,beq] = getConstraintOptionsLHFBA(model,lhfba_options)
    %Initialize
    Ncells = model.sizeCells;
    Nrxns = model.sizeYrxn + model.sizeXrxn + model.sizePrxn;
    Noptions = length(lhfba_options);
    Aineq = zeros(Ncells*Noptions,Ncells*(Nrxns+1));
    bineq = zeros(Ncells*Noptions,1);
    Aeq = zeros(Ncells*Noptions,Ncells*(Nrxns+1));
    beq = zeros(Ncells*Noptions,1);
    [ndf_idx,vxndf_idx] = getVariableIdxLHFBA(model);

    %Iterate over options
    Neqs = 0;
    Nineqs = 0;
    for option_idx = 1:Noptions
        option_entry = lhfba_options{option_idx};
        rxn_name = option_entry{1};
        bound = option_entry{2};
        value = option_entry{3};
        if length(value) == 1
            value = repmat(value,Ncells,1);
        end
        rxn_idx = strcmp(rxn_name,model.rxns);
        
        switch bound
            case 'l'
                %Iterate over cells
                for cell_idx = 1:Ncells
                    Aineq(Nineqs+cell_idx,ndf_idx(cell_idx)) = value(cell_idx)*model.xcells(cell_idx);
                    Aineq(Nineqs+cell_idx,vxndf_idx(rxn_idx,cell_idx)) = -1;
                end
                Nineqs = Nineqs + Ncells;
            case 'u'
                %Iterate over cells
                for cell_idx = 1:Ncells
                    Aineq(Nineqs+cell_idx,ndf_idx(cell_idx)) = -value(cell_idx)*model.xcells(cell_idx);
                    Aineq(Nineqs+cell_idx,vxndf_idx(rxn_idx,cell_idx)) = 1;
                end
                Nineqs = Nineqs + Ncells;
            otherwise
                %Iterate over cells
                for cell_idx = 1:Ncells
                    Aeq(Neqs+cell_idx,ndf_idx(cell_idx)) = value(cell_idx)*model.xcells(cell_idx);
                    Aeq(Neqs+cell_idx,vxndf_idx(rxn_idx,cell_idx)) = -1;
                end
                Neqs = Neqs + Ncells;
        end
    end
    
    %Remove unused equality/inequality space (done for matrix memory
    %preallocation)
    Aineq = Aineq(1:Nineqs,:);
    bineq = bineq(1:Nineqs,1);
    Aeq = Aeq(1:Neqs,:);
    beq = beq(1:Neqs,1);
end