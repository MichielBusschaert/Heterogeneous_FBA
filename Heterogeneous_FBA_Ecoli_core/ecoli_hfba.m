%% Clear
%clear all;
close all;
clc;

%% Generate HFBA model
%Load E. Coli model
fba_model = load('ecoli_fba_model.mat').model;

%Conversion from cell size to cell mass (assume constant diameter of 1 μm)
ecoli_diameter = 1;
ecoli_density = 1.105; %pgCDW/μm^3
x_ecoli = @(L) ecoli_density*pi*ecoli_diameter^2/4*L; %L: μm; x_ecoli: pgCDW
L_ecoli = @(x) x/ecoli_density/pi/ecoli_diameter^2*4; %x: pgCDW; L_ecoli: μm
L_max = 10 ; %μm

%Define heterogeneous grid
xmin = 0;
xmax = x_ecoli(L_max);
Ncells = 20;

%Division rate - constant
n_div = 12;
h_div = 5.65;
k_div = .155*60;
div_gamma_ct = @(x) k_div*(L_ecoli(x).^n_div)./(h_div.^n_div+L_ecoli(x).^n_div);

%Division rate - growth rate dependent
xdiv_min = .3;
xdiv_scale = 40;
% div_gamma_mu = @(x) max(xdiv_scale*(x-xdiv_min),0);
div_gamma_mu = @(x) zeros(size(x));

%Division kernel
bavg = 2;
a_inf = 1e5;
div_beta = @(x,y) bavg*sqrt(a_inf/pi)*exp(-(x-y/2).^2*a_inf);

%Coordinate transform
div_beta = @(z,w) div_beta(z/xmax,w/xmax)/xmax;

%Create HFBA model
hfba_model = generate_HFBA_model(fba_model,xmin,xmax,Ncells,div_gamma_ct,div_gamma_mu,div_beta);

%HFBA options
hfba_options = {{'ATPM','l',8.39};
    {'EX_glc(e)','l',-10};
    {'EX_ac(e)','l',0};
    {'EX_acald(e)','l',0};
    {'EX_akg(e)','l',0};
    {'EX_etoh(e)','l',0}
    {'EX_glu-L(e)','l',0};
    {'EX_pyr(e)','l',0};
    {'EX_lac-D(e)','l',0}};
% hfba_options = {{'ATPM','l',8.39};
%     {'EX_glc(e)','l',-10}};
hfba_model.uptakeLim(strcmp('o2[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = -12;
hfba_model.uptakeLim(strcmp('ac[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('acald[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('akg[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('etoh[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('glu-L[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('pyr[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;
hfba_model.uptakeLim(strcmp('lac-D[EXT]',hfba_model.mets(1:hfba_model.sizeYmet))) = 0;

% {'EX_ac(e)','l',0};
% {'EX_acald(e)','l',0};
% {'EX_akg(e)','l',0};
% {'EX_etoh(e)','l',0}
% {'EX_glu-L(e)','l',0};
% {'EX_pyr(e)','l',0};
% {'EX_lac-D(e)','l',0}}

%FBA options
fba_options = {{'ATPM','l',8.39};
    {'EX_glc(e)','l',-10};
    {'EX_o2(e)','l',-12}};

%% Solve regular FBA
v_fba = FBA_function(fba_model,fba_options,'cplex');
mu_fba = fba_model.c'*v_fba;
fba_table = table(hfba_model.rxns,rxn_equation(hfba_model.rxn_expr,hfba_model.rev),v_fba,'VariableNames',{'Rxn','Expr','Flux'});

%% Solve HFBA model
B0 = 1;
[ndf_sol,v_sol,mu_avg,sol_object] = LHFBA_distribution(hfba_model,mu_fba,B0,hfba_options);

% %% Visualize
% %NDF
% figure_ndf = figure();
% clf;
% plot(L_ecoli(hfba_model.xcells),ndf_sol,'b-','LineWidth',3);
% xlabel('Cell size (μm)','FontSize',14);
% ylabel('NDF (#/m^3/gCDW)','FontSize',14);
% title('Number Density Function','FontSize',16);
% 
% %Biomass
% figure_bio = figure();
% clf;
% plot(L_ecoli(hfba_model.xcells),hfba_model.xcells.*ndf_sol.*(hfba_model.ub-hfba_model.lb),'b-','LineWidth',3);
% xlabel('Cell size (μm)','FontSize',14);
% ylabel('Biomass (gCDW/m^3)','FontSize',14);
% title('Total biomass','FontSize',16);
% 
% %Growth rate
% figure_mu = figure();
% clf;
% hold on;
% plot(L_ecoli(hfba_model.xcells),v_sol'*hfba_model.c./B0,'b-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),repmat(mu_fba,hfba_model.sizeCells,1),'r--','LineWidth',2);
% plot(L_ecoli(hfba_model.xcells),repmat(mu_avg,hfba_model.sizeCells,1),'g--','LineWidth',2);
% xlabel('Cell size (μm)','FontSize',14);
% ylabel('\mu (1/h)','FontSize',14);
% title('Growth rate','FontSize',16);
% 
% %Glucose uptake
% figure_glc = figure();
% clf;
% hold on;
% plot(L_ecoli(hfba_model.xcells),v_sol(strcmp(hfba_model.rxns,'EX_glc(e)'),:),'b-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),repmat(v_fba(strcmp(hfba_model.rxns,'EX_glc(e)')),hfba_model.sizeCells,1),'r--','LineWidth',2);
% xlabel('Cell size (μm)','FontSize',14);
% ylabel('C6 Uptake (1/h)','FontSize',14);
% title('C6 uptake','FontSize',16);
% 
% %Oxygen uptake
% figure_O2 = figure();
% clf;
% hold on;
% plot(L_ecoli(hfba_model.xcells),v_sol(strcmp(hfba_model.rxns,'EX_glc(e)'),:),'b-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),repmat(v_fba(strcmp(hfba_model.rxns,'EX_glc(e)')),hfba_model.sizeCells,1),'r--','LineWidth',2);
% xlabel('Cell size (μm)','FontSize',14);
% ylabel('O2 Uptake (1/h)','FontSize',14);
% title('O2 uptake','FontSize',16);

%% Flux table
hfba_table = table(hfba_model.rxns,rxn_equation(hfba_model.rxn_expr,hfba_model.rev),v_sol,'VariableNames',{'Rxn','Expr','Flux per size range'});

%% Variability analysis
%Perform variability analysis
[ndf_lb,v_lb,ndf_ub,v_ub] = LHFBA_variability(hfba_model,mu_avg,B0,hfba_options,'cplex');

%NDF - Variability
figure_ndf_var = figure();
clf;
hold on;
plot(L_ecoli(hfba_model.xcells),ndf_sol*1e12,'b-','LineWidth',3);
fill(L_ecoli([hfba_model.xcells;flip(hfba_model.xcells)]),[ndf_lb;flip(ndf_ub)]*1e12,'b','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
%plot(L_ecoli(hfba_model.xcells),ndf_lb,'r-','LineWidth',3);
%plot(L_ecoli(hfba_model.xcells),ndf_ub,'r-','LineWidth',3);
xlabel('Cell length (μm)','FontSize',14);
ylabel('NDF (# gCDW^{-1})','FontSize',14);
%title('Number Density Function','FontSize',16);
grid on;

%Biomass - Variability
figure_bio_var = figure();
clf;
hold on;
plot(L_ecoli(hfba_model.xcells),hfba_model.xcells.*ndf_sol,'b-','LineWidth',3);
fill(L_ecoli([hfba_model.xcells;flip(hfba_model.xcells)]),[hfba_model.xcells.*ndf_lb;flip(hfba_model.xcells.*ndf_ub)],'b','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
% plot(L_ecoli(hfba_model.xcells),hfba_model.xcells.*ndf_lb.*(hfba_model.ub-hfba_model.lb),'r-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),hfba_model.xcells.*ndf_ub.*(hfba_model.ub-hfba_model.lb),'r-','LineWidth',3);
xlabel('Cell length (μm)','FontSize',14);
ylabel('Biomass density (gCDW pgCDW^{-1})','FontSize',14);
%title('Total biomass','FontSize',16);
grid on;

%Growth rate - Variability
figure_mu_var = figure();
clf;
hold on;
plot(L_ecoli(hfba_model.xcells),v_sol'*hfba_model.c./B0,'b-','LineWidth',3);
plot(L_ecoli(hfba_model.xcells),repmat(mu_avg,hfba_model.sizeCells,1),'r--','LineWidth',3);
mu_lb = v_lb'*hfba_model.c./B0;
mu_ub = v_ub'*hfba_model.c./B0;
vis_idces = ~(isnan(mu_lb)|isinf(mu_lb)|isnan(mu_ub)|isinf(mu_ub));
fill(L_ecoli([hfba_model.xcells(vis_idces);flip(hfba_model.xcells(vis_idces))]),[mu_lb(vis_idces);flip(mu_ub(vis_idces))],'b','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
% plot(hfba_model.xcells,mu_lb,'r-','LineWidth',3);
% plot(hfba_model.xcells,mu_ub,'r-','LineWidth',3);
xlabel('Cell length (μm)','FontSize',14);
ylabel('Growth rate \mu (h^{-1})','FontSize',14);
%title('Growth rate','FontSize',16);
legend({'HFBA','FBA'},'FontSize',12);
ylim([0,.8]);
grid on;

%Glucose uptake - Variability
figure_glc = figure();
clf;
hold on;
plot(L_ecoli(hfba_model.xcells),v_sol(strcmp(hfba_model.rxns,'EX_glc(e)'),:),'b-','LineWidth',3);
plot(L_ecoli(hfba_model.xcells),repmat(v_fba(strcmp(hfba_model.rxns,'EX_glc(e)')),hfba_model.sizeCells,1),'g--','LineWidth',3);
v_lb_glc = v_lb(strcmp(hfba_model.rxns,'EX_glc(e)'),:)';
v_ub_glc = v_lb(strcmp(hfba_model.rxns,'EX_glc(e)'),:)';
v_lb_glc(isnan(v_lb_glc)|isinf(v_lb_glc)) = 0;
v_ub_glc(isnan(v_ub_glc)|isinf(v_ub_glc)) = 0;
fill(L_ecoli([hfba_model.xcells;flip(hfba_model.xcells)]),[v_lb_glc;flip(v_ub_glc)],'b','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
% plot(L_ecoli(hfba_model.xcells),v_lb_glc,'r-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),v_ub_glc,'r-','LineWidth',3);
xlabel('Cell size (μm)','FontSize',14);
ylabel('C6 Uptake (1/h)','FontSize',14);
%title('C6 uptake','FontSize',16);
legend({'HFBA','FBA'},'FontSize',12);
grid on;

%Oxygen uptake - Variability
figure_O2 = figure();
clf;
hold on;
plot(L_ecoli(hfba_model.xcells),v_sol(strcmp(hfba_model.rxns,'EX_o2(e)'),:),'b-','LineWidth',3);
plot(L_ecoli(hfba_model.xcells),repmat(v_fba(strcmp(hfba_model.rxns,'EX_o2(e)')),hfba_model.sizeCells,1),'r--','LineWidth',3);
v_lb_o2 = v_lb(strcmp(hfba_model.rxns,'EX_o2(e)'),:)';
v_ub_o2 = v_lb(strcmp(hfba_model.rxns,'EX_o2(e)'),:)';
v_lb_o2(isnan(v_lb_o2)|isinf(v_lb_o2)) = 0;
v_ub_o2(isnan(v_ub_o2)|isinf(v_ub_o2)) = 0;
fill(L_ecoli([hfba_model.xcells;flip(hfba_model.xcells)]),[v_lb_o2;flip(v_ub_o2)],'b','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
% plot(L_ecoli(hfba_model.xcells),v_lb_o2,'r-','LineWidth',3);
% plot(L_ecoli(hfba_model.xcells),v_ub_o2,'r-','LineWidth',3);
xlabel('Cell size (μm)','FontSize',14);
ylabel('O2 Uptake (1/h)','FontSize',14);
%title('O2 uptake','FontSize',16);
legend({'HFBA','FBA'},'FontSize',12);
grid on;

%Table
flux_table_lb = table(hfba_model.rxns,rxn_equation(hfba_model.rxn_expr,hfba_model.rev),v_lb,'VariableNames',{'Rxn','Expr','Flux per size range'});
flux_table_ub = table(hfba_model.rxns,rxn_equation(hfba_model.rxn_expr,hfba_model.rev),v_ub,'VariableNames',{'Rxn','Expr','Flux per size range'});