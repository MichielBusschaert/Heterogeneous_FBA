%% Load e. coli model
clear all;
clc;
load('ecoli_core_model.mat');

%% Reorder metabolites
met_ext_id = contains(model.mets,'[e]');
met_met_id = contains(model.mets,'[c]');

mets_new = [model.mets(met_ext_id);model.mets(met_met_id)];
S_new = [model.S(met_ext_id,:);model.S(met_met_id,:)];
metCharge_new = [model.metCharge(met_ext_id);model.metCharge(met_met_id)];
metNames_new = [model.metNames(met_ext_id);model.metNames(met_met_id)];
metFormulas_new = [model.metFormulas(met_ext_id);model.metFormulas(met_met_id)];

%% Reorder reactions
rxn_ext_id = contains(model.rxns,'(e)');
rxn_bio_id = contains(model.rxns,'Biomass');
rxn_met_id = ~(rxn_ext_id|rxn_bio_id);

rxns_new = [model.rxns(rxn_ext_id);model.rxns(rxn_met_id);'biomass_production'];
S_new = [S_new(:,rxn_ext_id),S_new(:,rxn_met_id),S_new(:,rxn_bio_id)];
rev_new = [model.rev(rxn_ext_id);model.rev(rxn_met_id);model.rev(rxn_bio_id)];
c_new = [model.c(rxn_ext_id);model.c(rxn_met_id);model.c(rxn_bio_id)];
subSystems_new = [model.subSystems(rxn_ext_id);model.subSystems(rxn_met_id);model.subSystems(rxn_bio_id)];
rxnNames_new = [model.rxnNames(rxn_ext_id);model.rxnNames(rxn_met_id);model.rxnNames(rxn_bio_id)];

uptakeLim = zeros(sum(rxn_ext_id),1);
uptakeLim(model.lb(rxn_ext_id)~=0) = -Inf;

%% Add biomass metabolite
mets_new = [mets_new;'biomass'];
S_new = [S_new;c_new'];
metCharge_new = [metCharge_new;0];
metNames_new = [metNames_new;'Biomass'];
metFormulas_new = [metFormulas_new;{''}];

%% Add extracellular metabolites to exchange reactions
ex_mets_id = contains(mets_new,'[e]');
mets_new = [replace(mets_new(ex_mets_id),'[e]','[EXT]');mets_new];
metCharge_new = [metCharge_new(ex_mets_id);metCharge_new];
metNames_new = [append('EXT. ',metNames_new(ex_mets_id));metNames_new];
metFormulas_new = [metFormulas_new(ex_mets_id);metFormulas_new];

%% Update stoichiometric matrix
%Add ydot part
Sy = zeros(sum(ex_mets_id),size(S_new,2));
Sy(:,1:sum(ex_mets_id)) = -S_new(1:sum(ex_mets_id),1:sum(ex_mets_id));
S_new = [Sy;S_new];

%% Add mets info
model.sizeYmet = sum(met_ext_id);
model.sizeXmet = sum(met_ext_id)+sum(met_met_id);
model.sizePmet = 1;
model.sizeYrxn = sum(rxn_ext_id);
model.sizeXrxn = sum(rxn_met_id);
model.sizePrxn = 1;

%% Remove fields
model = rmfield(model,{'lb','ub','rules','genes','rxnGeneMat','grRules','confidenceScores','rxnReferences','rxnECNumbers','rxnNotes','metChEBIID','metKeggID','metPubChemID','metInchiString','b'});

%% Reassignment
model.rxns = rxns_new;
model.mets = mets_new;
model.S = S_new;
model.rev = rev_new;
model.c = c_new;
model.metCharge = metCharge_new;
model.subSystems = subSystems_new;
model.rxnNames = rxnNames_new;
model.metNames = metNames_new;
model.metFormulas = metFormulas_new;
model.uptakeLim = uptakeLim;

%% Add reaction expressions
model = createRxnExpression(model);

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