function [Aineq,bineq,Aeq,beq,lb,ub] = formulateConstraintsLHFBA(model,mu,B0,lhfba_options)      
    %Heterogeneity constraint
    [Aeq1,beq1] = getConstraintHeterogeneityLHFBA(model,mu);
    
    %Heterogeneity constraint
    [Aeq2,beq2] = getConstraintBoundaryConditionLHFBA(model);
    
    %Metabolism constraint
    [Aeq3,beq3] = getConstraintMetabolismLHFBA(model);
    
    %Exchange constraint
    [Aineq1,bineq1] = getConstraintExchangeLHFBA(model);
    
    %Biomass production constraint
    [Aeq4,beq4] = getConstraintBiomassProductionLHFBA(model,mu,B0);
    
    %Total biomass constraint
    [Aeq5,beq5] = getConstraintTotalBiomassLHFBA(model,B0);
    
    %Options constraint
    [Aineq2,bineq2,Aeq6,beq6] = getConstraintOptionsLHFBA(model,lhfba_options);
    
    %Get bounds
    [lb,ub] = getBoundsLHFBA(model);
    
    %Group constraints
%     Aineq = [Aineq1;Aineq2];
%     bineq = [bineq1;bineq2];
%     Aeq = [Aeq1;Aeq2;Aeq3;Aeq4;Aeq5;Aeq6];
%     beq = [beq1;beq2;beq3;beq4;beq5;beq6];
    Aineq = [Aeq1;-Aeq4;Aineq1;Aineq2];
    bineq = [beq1;-beq4;bineq1;bineq2];
    Aeq = [Aeq2;Aeq3;Aeq5;Aeq6];
    beq = [beq2;beq3;beq5;beq6];
end