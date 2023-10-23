function [obj,Rfirst,options] = initReach_inputDependence(obj,Rinit,options)
% initReach_inputDependence - computes the continuous reachable continuous 
%    for the first time step when the constant input is parameterized and
%    correlated to the parameters of the system
%
% Syntax:
%    [obj,Rfirst,options] = initReach_inputDependence(obj,Rinit,options)
%
% Inputs:
%    obj - linParamSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linParamSys object
%    Rfirst - first reachable set 
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
% Written:       01-June-2011
% Last update:   15-February-2021 (MW, rename: intermediateTerms)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% store taylor terms and time step as object properties
obj.stepSize = options.timeStep;
obj.taylorTerms = options.taylorTerms;

% compute mapping matrix
obj = mappingMatrix(obj,options);
% compute high order mapping matrix
obj = highOrderMappingMatrix(obj,options.intermediateTerms);
% compute time interval error (tie)
obj = tie(obj);
% compute reachable set due to input
obj = inputSolution(obj,options);

%compute reachable set of first time interval
%first time step homogeneous solution
Rhom_tp = dependentHomSol(obj, Rinit, options.Uconst);

%time interval solution
inputCorr = obj.inputF*obj.B*zonotope(options.uTrans + center(options.Uconst));
Rhom = enclose(Rinit,Rhom_tp) + obj.F*Rinit + inputCorr;

%total solution
Rtotal = Rhom + obj.RV;
Rtotal_tp = Rhom_tp + obj.RV;

%write results to reachable set struct Rfirst
Rfirst.tp = reduce(Rtotal_tp,options.reductionTechnique,options.zonotopeOrder);
Rfirst.ti = reduce(Rtotal,options.reductionTechnique,options.zonotopeOrder);

% ------------------------------ END OF CODE ------------------------------
