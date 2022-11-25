function [Rfirst,options] = initReach(obj, Rinit, options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:  
%    [obj,Rfirst] = initReach(obj,Rinit,options)
%
% Inputs:
%    obj - linParamSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linParamSys object
%    Rfirst - first reachable set 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-August-2010
% Last update:  16-May-2011
%               19-February-2012
%               12-August-2016
%               19-May-2020 (MW, error handling for exploding sets)
% Last revision:---

%------------- BEGIN CODE --------------

% store taylor terms and time step as object properties
obj.stepSize = options.timeStep;
obj.taylorTerms = options.taylorTerms;

% compute mapping matrix
mappingMatrix(obj,options);
% compute time interval error (tie)
tie(obj);
% compute reachable set due to input
inputSolution(obj,options);
%change the time step size
obj.stepSize=options.timeStep;

%compute reachable set of first time interval
%first time step homogeneous solution
Rhom_tp = obj.mappingMatrixSet.zono*Rinit + obj.mappingMatrixSet.int*Rinit;
if isa(Rinit,'zonoBundle')
    Rhom = enclose(Rinit,Rhom_tp+obj.Rtrans) ...
    + obj.F*Rinit.Z{1} + obj.inputCorr + (-1*obj.Rtrans);
else
    Rhom = enclose(Rinit,Rhom_tp+obj.Rtrans) ...
    + obj.F*Rinit + obj.inputCorr + (-1*obj.Rtrans);
end

% total solution
Rtotal = Rhom + obj.Rinput;
Rfirst.ti = reduce(Rtotal,options.reductionTechnique,options.zonotopeOrder);

if options.compTimePoint
    Rtotal_tp = Rhom_tp + obj.Rinput;
    Rfirst.tp = reduce(Rtotal_tp,options.reductionTechnique,options.zonotopeOrder);
else
    Rfirst.tp = [];
end

% check for explosion
if isa(Rinit,'zonotope')
    temp = rad(interval(Rinit));
    if all(temp > eps) && options.compTimePoint && ...
       max(rad(interval(Rfirst.tp)) ./ rad(interval(Rinit))) > 1e10 % arbitary value
        throw(printExplosionError());
    end
end


%------------- END OF CODE --------------