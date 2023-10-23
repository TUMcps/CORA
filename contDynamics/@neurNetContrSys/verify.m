function [res, R, simRes] = verify(obj, spec, params, options, evParams, varargin)
% verify - tries to verify the specification with the given params/options
%    1. simulation random rusn
%    2. check violations in simulations
%       2.1 if violated: set initial set to starting point
%       2.2 no violation: keep initial set
%    3. compute reachable set
%
% Syntax:
%    [res, R, simRes] = verify(obj, spec, params, options, evParams)
%    res = verify(obj, spec, params, options, evParams)
%
% Inputs:
%    obj - neurNetContrSys object
%    spec - object of class specification
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
%    evParams - parameters for neural network evaluation
%    verbose - verbose output
%
% Outputs:
%    res - 'VIOLATED', 'VERIFIED', 'UNKNOWN'
%    R - object of class reachSet storing the computed reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys

% Authors:       Tobias Ladner
% Written:       22-March-2023
% Last update:   20-July-2023 (TL, check early stop)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 5
    throw(CORAerror("CORA:notEnoughInputArgs", 5))
elseif nargin > 6
    throw(CORAerror('CORA:tooManyInputArgs', 6))
end
verbose = setDefaultValues({false}, varargin);

inputArgsCheck({ ...
    {obj, 'att', 'neurNetContrSys'}; ...
    {spec, 'att', 'specification'}; ...
    {params, 'att', 'struct'}; ...
    {options, 'att', 'struct'}; ...
    {evParams, 'att', 'struct'}; ...
    {verbose, 'att', 'logical', 'scalar'}; ...
})

% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(obj, params);
tSim = toc;
if verbose
    disp(['Time to compute random simulations: ', num2str(tSim)]);
end

% Check Violation ---------------------------------------------------------

tic
[isNotVio, ~, indObj] = check(spec, simRes);
isVio = ~isNotVio; % spec is safeSet
tVio = toc;
if verbose
    disp(['Time to check violation in simulations: ', num2str(tVio)]);
end

if isVio
    % only continue with violating run
    simRes = simResult( ...
        simRes(indObj{1}).x(indObj{2}), ...
        simRes(indObj{1}).t(indObj{2}));
    params.R0 = polyZonotope(simRes.x{1}(1, :)');
    
end

% Reachability Analysis ---------------------------------------------------

tic
R = reach(obj, params, options, evParams, spec);
tComp = toc;
if verbose
    disp(['Time to compute reachable set: ', num2str(tComp)]);
end

RtFinal = R(end).timePoint.time{end};
if RtFinal ~= params.tFinal
    fprintf("Reachability analysis stopped early at t=%g.\n",RtFinal)
    res = 'UNKNOWN';
    return
end

% Verification ------------------------------------------------------------

tic
isVeri = check(spec, R);
tVeri = toc;
if verbose
    disp(['Time to check verification: ', num2str(tVeri)]);
end

% Finish ------------------------------------------------------------------

if isVio && ~isVeri
    res = 'VIOLATED';
elseif isVeri
    res = 'VERIFIED';
else
    res = 'UNKNOWN';
end

if verbose
    tTotal = tSim+tVio+tComp+tVeri;
    disp(['Total Time: ', num2str(tTotal)]);
end

end

% ------------------------------ END OF CODE ------------------------------
