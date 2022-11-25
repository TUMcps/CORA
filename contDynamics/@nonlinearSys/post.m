function [Rnext,options] = post(obj,R,options)
% post - computes the reachable continuous set for one time step of a
% nonlinear system by overapproximative linearization
%
% Syntax:  
%    [Rnext] = post(obj,R,options)
%
% Inputs:
%    obj - nonlinearSys object
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
% Written:      03-January-2008
% Last update:  29-June-2009
%               10-August-2016
%               19-November-2017
% Last revision: ---

%------------- BEGIN CODE --------------

% potentially restructure the polynomial zonotope
if strcmp(options.alg,'poly') && isa(R.tp{1}.set,'polyZonotope') && ...
   isfield(options,'polyZono') && ~isinf(options.polyZono.maxPolyZonoRatio)

    temp = options.polyZono;

    for i = 1:length(R.tp)
        
        % compute ratio of dependent to independent part 
        ratio = approxVolumeRatio(R.tp{i}.set,temp.volApproxMethod);

        % restructure the polynomial zonotope
        if ratio > temp.maxPolyZonoRatio
           R.tp{i}.set = restructure(R.tp{i}.set, ...
                      temp.restructureTechnique,temp.maxDepGenOrder);
        end
    end
end

%In contrast to the linear system: the nonlinear system has to be constantly
%initialized due to the linearization procedure
[Rnext,options] = initReach(obj,R.tp,options);

%reduce zonotopes
for i=1:length(Rnext.tp)
    if ~isempty(Rnext.tp{i}.set)
        Rnext.tp{i}.set=reduce(Rnext.tp{i}.set,options.reductionTechnique,options.zonotopeOrder);
        Rnext.ti{i}=reduce(Rnext.ti{i},options.reductionTechnique,options.zonotopeOrder);
    end
end

%delete redundant reachable sets
Rnext = deleteRedundantSets(Rnext,R,options);

%------------- END OF CODE --------------