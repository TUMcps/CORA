function [Rnext,options] = post(obj,R,options)
% post - computes the reachable continuous set for one time step of a
% linear interval system
%
% Syntax:  
%    [Rnext] = post(obj,R,options)
%
% Inputs:
%    obj - linIntSys object
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

% Author:       Matthias Althoff
% Written:      15-May-2007 
% Last update:  07-January-2009
%               22-June-2009
%               29-June-2009
%               06-August-2010
%               08-August-2016
%               09-February-2016
% Last revision: ---

%------------- BEGIN CODE --------------

% rename for readabilitry
M1 = obj.mappingMatrixSet.zono;
M2 = obj.mappingMatrixSet.int;

% next set (no input)
R.ti = M1*R.ti + M2*R.ti; 

% bloating due to input
R.ti = R.ti + obj.Rinput;

% reduce zonotope
Rred=reduce(R.ti,options.reductionTechnique,options.zonotopeOrder);

% write results to reachable set struct Rnext
Rnext.ti=Rred;

if isfield(options,'includeTimePointSolution') && (options.includeTimePointSolution==1)
    % next set (no input)
    R.tp = M1*R.tp + M2*R.tp;

    % bloating due to input
    R.tp = R.tp + obj.Rinput;

    % reduce zonotope
    Rred_tp=reduce(R.tp,options.reductionTechnique',options.zonotopeOrder);

    % write results to reachable set struct Rnext
    Rnext.tp=Rred_tp;
else
    Rnext.tp = [];
end


%------------- END OF CODE --------------