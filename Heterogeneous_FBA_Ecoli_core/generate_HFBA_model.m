function [hfba_model] = generate_HFBA_model(fba_model,xmin,xmax,Ncells,division_gamma_ct,division_gamma_mu,division_beta)
    %Copy model
    hfba_model = fba_model;
    
    %Define grid points of the heterogeneous system
    grid_points = linspace(xmin,xmax,Ncells+1)';
    hfba_model.lb = grid_points(1:Ncells);
    hfba_model.ub = grid_points(2:Ncells+1);
    hfba_model.xcells = .5*(hfba_model.lb+hfba_model.ub);
    hfba_model.sizeCells = Ncells;
    
    %Assign division parameters
    hfba_model.division_gamma_ct = division_gamma_ct;
    hfba_model.division_gamma_mu = division_gamma_mu;
    hfba_model.division_beta = division_beta;
    hfba_model.division_beta_mat = division_kernel_average(hfba_model);
end