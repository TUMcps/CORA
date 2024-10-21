function [Rnext,options] = post(nlnsysDA,R,params,options)
% post - computes the reachable continuous set for one time step of a
%    nonlinear differential-algebraic system by over-approximative
%    abstraction
%
% Syntax:
%    [Rnext,options] = post(nlnsysDA,R,params,options)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%    R - reachable set of the previous time step
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
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
[Rnext,options] = initReach(nlnsysDA,R,params,options);

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
