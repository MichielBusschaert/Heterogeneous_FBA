function [model] = createRxnExpression(model)
    %Initialize
    Nrxns = length(model.rxns);
    rxn_expr = cell(Nrxns,1);
    
    %Iterate over reactions
    for rxn_idx = 1:Nrxns
        %Rxn cell
        stoich_id = model.S(:,rxn_idx)~=0;
        rxn_cell = cell(1,2*(sum(stoich_id)));
        mets_list = model.mets(stoich_id);
        coeff_list = model.S(stoich_id,rxn_idx);
        for met_idx = 1:sum(stoich_id)
            rxn_cell{2*met_idx-1} = mets_list{met_idx};
            rxn_cell{2*met_idx} = coeff_list(met_idx);
        end
        rxn_expr{rxn_idx} = rxn_cell;
    end
    model.rxn_expr = rxn_expr;
end