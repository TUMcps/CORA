function [Rnext,options] = post_Krylov(obj,options)
% post_Krylov - computes the reachable continuous set for one time step in 
% the Krylov subspace
%
% Syntax:  
%    [Rnext,options] = post_Krylov(obj,R,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-December-2016
% Last update:  03-March-2017
%               21-August-2017
%               02-November-2018
% Last revision:---

%------------- BEGIN CODE --------------

% retrieve reachable set
Rpar_proj = options.Rpar_proj;
R_tp_prev = options.Rhom_tp_proj;

% Multiply previous reachable set with exponential matrix
% retrieve initial reachable sets
%R_hom_0 = options.Rhom_0;
RV_0 = options.RV_0;
% multiplications with exponential matrix
[R_hom_tp_proj,R_Krylov_proj] = exponential_Krylov_projected(obj,options.R0,options,1);
[RV_proj,R_V_Krylov_proj] = exponential_Krylov_projected(obj,RV_0,options,0);

% other results
R_tie_proj = options.R_tie_proj;
inputCorr_proj = options.inputCorr_proj;
Rtrans_proj = options.Rtrans_proj;
%Rpar_proj = Rpar_proj + interval(RV_proj);
Rpar_proj = Rpar_proj + RV_proj;
Rpar_proj = reduce(Rpar_proj,options.reductionTechnique,options.zonotopeOrder);

% next time step homogeneous solution
% Comment: Rtrans only considers the contant input for one time interval; 
% other time intervals are considered in input solution
Rhom_tp_proj = R_hom_tp_proj + R_Krylov_proj + Rtrans_proj; 
R_tp_proj = Rhom_tp_proj + zonotope(Rpar_proj);
R_ti_proj = enclose(R_tp_prev,R_tp_proj) + R_tie_proj + zonotope(Rpar_proj) + inputCorr_proj;

% order reduction
Rnext.ti = reduce(R_ti_proj,options.reductionTechnique,options.zonotopeOrder);
Rnext.tp = reduce(R_tp_proj,options.reductionTechnique,options.zonotopeOrder);

% update options
options.Rhom_tp_proj = Rhom_tp_proj;
options.Raux_proj = RV_proj;
options.Rpar_proj = Rpar_proj;

%------------- END OF CODE --------------