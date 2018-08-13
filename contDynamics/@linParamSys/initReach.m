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
% Last revision:---

%------------- BEGIN CODE --------------

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
Rhom = enclose(Rinit,Rhom_tp+obj.Rtrans) ...
    + obj.F*Rinit + obj.inputCorr + (-1*obj.Rtrans);

%total solution
Rtotal = Rhom + obj.Rinput;
Rtotal_tp = Rhom_tp + obj.Rinput;

%write results to reachable set struct Rfirst
Rfirst.tp = reduce(Rtotal_tp,'girard',options.zonotopeOrder);
Rfirst.ti = reduce(Rtotal,'girard',options.zonotopeOrder);


%------------- END OF CODE --------------