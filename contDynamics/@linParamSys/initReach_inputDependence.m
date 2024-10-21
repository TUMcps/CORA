function [sys,Rfirst,options] = initReach_inputDependence(sys,Rinit,params,options)
% initReach_inputDependence - computes the continuous reachable continuous 
%    for the first time step when the constant input is parameterized and
%    correlated to the parameters of the system
%
% Syntax:
%    [sys,Rfirst,options] = initReach_inputDependence(sys,Rinit,params,options)
%
% Inputs:
%    sys - linParamSys object
%    Rinit - initial reachable set
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    sys - linParamSys object
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       01-June-2011
% Last update:   15-February-2021 (MW, rename: intermediateTerms)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% store taylor terms and time step as object properties
sys.stepSize = options.timeStep;
sys.taylorTerms = options.taylorTerms;

% compute mapping matrix
sys = mappingMatrix(sys,params,options);
% compute high order mapping matrix
sys = highOrderMappingMatrix(sys,options.intermediateTerms);
% compute time interval error (tie)
sys = tie(sys);
% compute reachable set due to input
sys = inputSolution(sys,params,options);

%compute reachable set of first time interval
%first time step homogeneous solution
Rhom_tp = dependentHomSol(sys, Rinit, params.Uconst);

%time interval solution
inputCorr = sys.inputF*sys.B*zonotope(params.uTrans + center(params.Uconst));
Rhom = enclose(Rinit,Rhom_tp) + sys.F*Rinit + inputCorr;

%total solution
Rtotal = Rhom + sys.RV;
Rtotal_tp = Rhom_tp + sys.RV;

%write results to reachable set struct Rfirst
Rfirst.tp = reduce(Rtotal_tp,options.reductionTechnique,options.zonotopeOrder);
Rfirst.ti = reduce(Rtotal,options.reductionTechnique,options.zonotopeOrder);

% ------------------------------ END OF CODE ------------------------------
