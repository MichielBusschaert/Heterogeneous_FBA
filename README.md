# Overview of this repository
This repository is associated to a manuscript [1]. In the manuscript, we introduced an extension to Flux Balance Analysis, in which we include cell size heterogeneity. We illustrated our method on a small-scale *E. coli* metabolic network (72 metabolites, 95 reactions), as described in [2]. This repository contains the Matlab code that was used to simulate the numerical values reported in Section IV of [1], and the Figures in [1]. The code in this repository can be used to recreate the results reported in [1] applied to the *E. coli* model, the same heterogeneous FBA method can be applied to a custom user model. See the section **HFBA model construction** in this text for further info on how to construct such model.

## HFBA model construction
The heterogeneous FBA model is represented as a `stuct` object in Matlab. The `struct` ***requires*** following fields
* `sizeYmet`: Number of extracellular metabolites $n_Y$
* `sizeXmet`: Number of intracellular metabolites $n_X$
* `sizePmet`: Number of biomass components $n_P$ (often equal to 1)
* `sizeYrxn`: Number of exchange reactions $n_{r,Y}$
* `sizeXrxn`: Number of intracellular reactions $n_{r,X}$
* `sizePrxn`: Number of biomass production reactions $n_{r,P}$
* `sizeCells`: Number of discretization intervals, i.e. $n_{disc} = L/\Delta x$
* `S`: Stoichiometric matrix, such that `S(i,j)` contains the reaction coefficient of metabolite `i` in reaction `j` (dimension: $`(n_Y+n_X+n_P)\times (n_{r,Y}+n_{r,X}+n_{r,P})`$ )
* `rev`: Vector with reaction reversibility data, such that `rev(i)` contains a `1` if the corresponding reaction is reversible, `0` if it is irreversible (length: $`(n_{r,Y}+n_{r,X}+n_{r,P})`$.
* `c`: Vector or matrix containing weigth coefficients to determine the specific growth rate. If vector, this is implies a constant biomass composition over all cell sizes (length: $`(n_{r,Y}+n_{r,X}+n_{r,P}`$). If matrix, the cell size dependency on biomass composition is taken into account (dimension: $`n_{disc}\times(n_{r,Y}+n_{r,X}+n_{r,P}`$)
* `uptakeLim`: The total uptake limit for extracellular components, function of **b** in (12e) of [1] (length: $n_Y$).
* `x`: Cell mass value in the center of a discretization interval (length: $n_{disc}$).
* `lb`: Lower bound of the cell mass interval (length: $n_{disc}$).
* `ub`: Upper bound of the cell mass interval, typicall `ub(i)` = `lb(i+1)` (length: $n_{disc}$).
* `division_gamma`: Matlab anonymous function describing division rate $\gamma(x)$ from (12b) in [1] (definition: `@(x) gamma(x)`).
* `division_beta`: Matlab anonymous function describing divison kernel $\beta(x,x')$ from (12b) in [1] (definition: `@(x,y) beta(x,y)`).

In addition, the `struct` can contain ***optional*** fields
* `mets`: Cell vector with names or identifiers for each metabolite (length: $`(n_Y+n_X+n_P)`$).
* `rxns`: Cell vector with names or identifiers for each reaction (length: $`(n_{r,Y}+n_{r,X}+n_{r,P})`$).
* `metNames`: Cell vector with full metabolite names (length: $`(n_Y+n_X+n_P)`$).
* `rxnNames`: Cell vector with full reaction names (length: $`(n_{r,Y}+n_{r,X}+n_{r,P})`$).

## Content of this repository
The main code to execute the heterogeneous FBA method on an existing model is contained in the LHFBA folder. Below, the function of each script is briefly explained. Note that to run the script, the user needs to have either ILOG IBM CPLEX installed, or the Matlab Optimization Toolbox (which contains the function `linprog`). The results in [1] were generated using CPLEX.

### LHFBA
Following files are included in the LHFBA folder, and are sufficient to obtain results given an HFBA model `struct`:
* `LHFBA_distribution.m`: This is the main function, which should be called by the user to solve the heterogeneous FBA program. An initial guess `mu_avg_init` should be supplied. In addition, `hfba_options` contains info on upper and lower reaction flux bounds $`v_{lb}(x)`$ and $`v_{ub}(x)`$ in (12h) in [1]. This option is defined as a cell array of cell array, where the cell array contains three entries, 1) the reaction index or the name of the reaction, 2) a character to denote lower, upper or equal bound ('l'/'u'/'e'), respectively, and 3) the value of the bound.
* `LHFBA_function.m`: Evaluates the program (12) in [1], and if feasible, returns the number density function and different fluxes, and if infeasible, returns empty vectors.
* `LHFBA_variability.m`: Performs a variability analysis; for each optimization variables, it calculates its upper and lower value for a fixed growth rate.
* `formulateConstraintsLHFBA.m`: Returns matrices `A`and `b` to impose the linear constraints, in the form of `A*x<=b` or `A*x=b`.
* `formulateObjectiveLHFBA.m`: Returns a vector with coefficients in the linear objective, `J = c'*x`. Since the feasibility is assessed, this objective could be generic.
* `getBoundsLHFBA.m`: Returns lower and upper bounds on the optimization variables, `lb <= x <= ub` (includes Equation (12i) in [1]).
* `getConstraintBiomassProductionLHFBA.m`: Returns the biomass production constraint (Equation (12f) in [1])
* `getConstraintBoundaryConditionLHFBA.m`: Returns the boundary condition on the PBE (Equation (12c) in [1])
* `getConstraintExchangeLHFBA.m`: Returns the extracellular mass balance constraint (Equation (12e) in [1])
* `getConstraintHeterogeneityLHFBA.m`: Returns the (discretized) population balance equation constraint (Equation (12b) in [1])
* `getConstraintMetabolismLHFBA.m`: Returns the intracellular mass balance constraint (Equation (12d) in [1])
* `getConstraintOptionsLHFBA.m`: Returns the flux constraints in terms of lower and upper bounds (Equation (12h) in [1])
* `getConstraintTotalBiomassLHFBA.m`: Returns the normalization constraint for the NDF (Equation (12g) in [1])
* `getVariableIdxLHFBA.m`: Helper function to easily access optimization variable indices.

### Other files
Following other files are included in this folder:
* The FBA directory: This directory contains a custom implementation to standard Flux Balance Analysis, which was used to compare results of HFBA to standard FBA.
* `convert_ecoli_model.m`: A helper file that specifically restructures the *E. coli* model.
* `createRxnExpression.m`: A helper file to construct character-based reaction expressions based on the stoichiometric matrix `S`.
* `ecoli_hfba.m`: The main script which constructs the HFBA model and solves the optimization program, as well as visualises the results.
* `generate_HFBA_model.m`: A helper function to convert a standard FBA model into a heterogeneous FBA model.

## *E. coli* model
Download *E. coli* model [2] from https://systemsbiology.ucsd.edu/Downloads/EcoliCore.

## How to run heterogeneous FBA
Given a heterogeneous FBA model as Matlab `struct`, simply call the function `LHFBA_distribution(model,mu_avg_init,biomass,hfba_options)`, where
* `model` is the HFBA model as Matlab `struct`, containing the fields given earlier,
* `mu_avg_init` is the initial guess for the average specific growth rate $`\bar{\mu}`$,
* `biomass` is the total biomass ($`B`$ in [1]),
* `hfba_options` is a cell array of 3-by-1 cell arrays containing info on flux bounds.
The output of this function is `ndf_sol`, `v_sol` and `mu_avg_sol`, where
* `ndf_sol` is a vector containing the NDF solution $`n(x)`$ (length: $`n_{disc}`$),
* `v_sol` is a matrix containing the reaction flux solutions $`v(x)`$ (dimension: $`(n_{r,Y}+n_{r,X}+n_{r,P}) \times n_{disc}`$,
* `mu_avg_sol` is the optimal average specific growth rate $`\bar{\mu}^*`$ (scalar)

## Info on discretization scheme and solution method
A discretization of \eqref{eq:hfba} is given below. A cell size range is defined from $0$ to a (sufficiently large) $L_\text{max}$ and is discretized into $N$ equally spaced intervals with size $\Delta x = L_\text{max}/N$. All cell sizes between the interval bounds $x_i \pm \Delta x/2$ are approximated with $x_i$ for $i\in \{1,2\ldots,N\}$. Similarly, the NDF $n(x)$ and the term $\mathbf{v}(x)xn(x)$ in the $i^\text{th}$ interval are represented as $n_i$ and $\mathbf{v}xn_i$. The PDE (12b) in [1] is integrated over each interval $[x_i-\Delta x/2,x_i+\Delta x/2]$. This results in evaluation of the derivatives to cell mass $x$ at the interval boundaries. Using a first order upwind scheme [3], one obtains:
```math
    \begin{align}
        \max_{n_i,\mathbf{v}xn_i} & \bar{\mu}, \\
        \text{s.t. } & \bar{\mu}n_i \leq - \frac{1}{\Delta x}\left(\mathbf{c}_i^\intercal\mathbf{v}xn_i-\mathbf{c}_{i-1}^\intercal\mathbf{v}xn_{i-1}\right) -\gamma(x_i)n_i + \sum_{j=i}^N \beta_{i,j}\gamma(x_j)n_j\Delta x, \\
        & \mathbf{c}_0^\intercal\mathbf{v}xn_0=\mathbf{c}_N^\intercal\mathbf{v}xn_N = 0, \\
        & S_{X}\cdot \mathbf{v}xn_i = \mathbf{0}, \\
        & \sum_{j=1}^N S_{Y}\cdot \mathbf{v}xn_j\Delta x \geq \mathbf{b}, \\
        & \sum_{j=1}^N \mathbf{c}_j^\intercal \mathbf{v}xn_j\Delta x \geq \bar{\mu}B, \\
        & \sum_{j=1}^N x_jn_j\Delta x = B, \\
        & \mathbf{v}_\text{lb}(x_i)x_in_i \leq \mathbf{v}xn_i \leq \mathbf{v}_\text{ub}(x_i)x_in_i, \\
        & 0 \leq n_i, \text{for all } i\in \{1,2,\ldots,N\}.
    \end{align}
```
In \eqref{seq:hfbad_b}, the breakage kernel $\beta(x,x')$, where $x\in$ interval $i$ and $x'\in$ interval $j$, is discretized as $\beta_{i,j}$ with
```math
\beta_{i,j} = \int_{x_i-\frac{\Delta x}{2}}^{\min\left(x_j,x_i+\frac{\Delta x}{2}\right)} \beta(\xi,x_j)d\xi.
```

# References

[1] M. Busschaert, F. H. Vermeire, and S. Waldherr, "Modeling Cell Size Distribution with Heterogeneous Flux Balance Analysis," *arXiv*, 2023. Available online: https://arxiv.org/abs/2304.06631.

[2] J. D. Orth, R. M. Fleming, and B. Ø. Palsson, “Reconstruction and use of microbial metabolic networks: the core *Escherichia coli* metabolic model as an educational guide,” *EcoSal plus*, vol. 4, no. 1, 2010.

[3] S. Motz, A. Mitrovi, and E.-D. Gilles, “Comparison of numerical methods for the simulation of dispersed phase systems,” *Chemical Engineering Science*, vol. 57, no. 20, pp. 4329–4344, 2002.
