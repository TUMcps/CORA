function [Rnext,options] = priv_post_Krylov_nonProjected(linsys,R,options)
% priv_post_Krylov_nonProjected - computes the reachable continuous set for one
%    time step in the Krylov subspace
%
% Syntax:
%    [Rnext,options] = priv_post_Krylov_nonProjected(linsys,R,options)
%
% Inputs:
%    linsys - linearSys object
%    R - reachable set of the previous time step
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

% Authors:       Matthias Althoff
% Written:       22-December-2016
% Last update:   03-March-2017
%                21-August-2017
%                02-November-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% retrieve last reachable set
Rpar = options.Rpar;

% Multiply previous reachable set with exponential matrix
% retrieve initial reachable sets
%R_hom_0 = options.Rhom_0;
RV_0 = options.RV_0;
% multiplications with exponential matrix
[R_hom_tp,R_Krylov] = priv_exponential_Krylov_precomputed(linsys,options.R0,options,1);
[RV,R_V_Krylov] = priv_exponential_Krylov_precomputed(linsys,RV_0,options,0);

% %COMMENT: old solution without precomputation
% % retrieve last reachable sets
% Rhom_tp = options.Rhom_tp;
% RV = options.Raux;
% % multiplications with exponential matrix
% [R_hom_tp,R_Krylov] = priv_exponential_Krylov(linsys,Rhom_tp,options);
% [RV,R_V_Krylov] = priv_exponential_Krylov(linsys,RV,options);


% other results
R_tie = options.R_tie;
inputCorr = options.inputCorr;
Rtrans = options.Rtrans;
Rpar = Rpar + interval(RV); % + interval(Rtrans) + (-center(Rtrans));

% next time step homogeneous solution
Rhom_tp = R_hom_tp + R_Krylov + Rtrans; %Rtrans not yet changed
R_tp = Rhom_tp + zonotope(Rpar);
R_ti = enclose(R.tp,R_tp) + R_tie + zonotope(Rpar) + inputCorr;

% order reduction
Rnext.ti = reduce(R_ti,options.reductionTechnique,options.zonotopeOrder);
Rnext.tp = reduce(R_tp,options.reductionTechnique,options.zonotopeOrder);

% update options
options.Rhom_tp = Rhom_tp;
options.Raux = RV;
options.Rpar = Rpar;

% ------------------------------ END OF CODE ------------------------------
