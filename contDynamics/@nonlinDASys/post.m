function [Rnext,options] = post(obj,R,options)
% post - computes the reachable continuous set for one time step of a
%    nonlinear differential-algebraic system by over-approximative
%    abstraction
%
% Syntax:
%    [Rnext,options] = post(obj,R,options)
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

% Authors:       Matthias Althoff
% Written:       22-November-2011
% Last update:   08-August-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%despite the linear system: the nonlinear system has to be constantly
%initialized due to the linearization procedure
[Rnext,options] = initReach(obj,R,options);

%reduce zonotopes
for i=1:length(Rnext.tp)
    % time-point solution
    Rnext.tp{i}.set = reduce(Rnext.tp{i}.set,...
        options.reductionTechnique,options.zonotopeOrder);
    % time-interval solution
    Rnext.ti{i} = reduce(Rnext.ti{i},...
        options.reductionTechnique,options.zonotopeOrder);
end

%delete redundant reachable sets
Rnext = deleteRedundantSets(Rnext,R,options);

% ------------------------------ END OF CODE ------------------------------
