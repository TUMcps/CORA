function [Rfirst,options] = initReach(sys,Rinit,params,options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:
%    [Rfirst,options] = initReach(sys,Rinit,params,options)
%
% Inputs:
%    sys - linParamSys object
%    Rinit - initial reachable set
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-August-2010
% Last update:   16-May-2011
%                19-February-2012
%                12-August-2016
%                19-May-2020 (MW, error handling for exploding sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% store taylor terms and time step as object properties
sys.stepSize = options.timeStep;
sys.taylorTerms = options.taylorTerms;

% compute mapping matrix
priv_mappingMatrix(sys,params,options);
% compute time interval error (tie)
priv_tie(sys);
% compute reachable set due to input
priv_inputSolution(sys,params,options);
%change the time step size
sys.stepSize=options.timeStep;

%compute reachable set of first time interval
%first time step homogeneous solution
Rhom_tp = sys.mappingMatrixSet.zono*Rinit + sys.mappingMatrixSet.int*Rinit;
if isa(Rinit,'zonoBundle')
    Rhom = enclose(Rinit,Rhom_tp+sys.Rtrans) ...
    + sys.F*Rinit.Z{1} + sys.inputCorr + (-1*sys.Rtrans);
else
    Rhom = enclose(Rinit,Rhom_tp+sys.Rtrans) ...
    + sys.F*Rinit + sys.inputCorr + (-1*sys.Rtrans);
end

% total solution
Rtotal = Rhom + sys.Rinput;
Rfirst.ti = reduce(Rtotal,options.reductionTechnique,options.zonotopeOrder);

if options.compTimePoint
    Rtotal_tp = Rhom_tp + sys.Rinput;
    Rfirst.tp = reduce(Rtotal_tp,options.reductionTechnique,options.zonotopeOrder);
else
    Rfirst.tp = [];
end

% check for explosion
if isa(Rinit,'zonotope')
    temp = rad(interval(Rinit));
    if all(temp > eps) && options.compTimePoint && ...
       max(rad(interval(Rfirst.tp)) ./ rad(interval(Rinit))) > 1e10 % arbitary value
        throw(CORAerror('CORA:reachSetExplosion'));
    end
end

% ------------------------------ END OF CODE ------------------------------
