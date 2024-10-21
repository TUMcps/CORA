function [R,res] = reach(linsys,params,varargin)
% reach - computes the reachable set for linear systems
%
% Syntax:
%    R = reach(linsys,params)
%    R = reach(linsys,params,options)
%    [R,res] = reach(linsys,params,options,spec)
%
% Inputs:
%    linsys - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true/false whether specifications are satisfied
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-October-2019
%                23-April-2020 (added params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
[options,spec] = setDefaultValues({struct('linAlg','adaptive'),[]},varargin);

% options preprocessing
[params,options] = validateOptions(linsys,params,options);

specLogic = [];
if ~isempty(spec)
    [spec,specLogic] = splitLogic(spec);
    if ~isempty(spec)
	    options.specification = spec;
    end
end

% hybrid systems: if invariant is empty set (used to model instant
% transitions), exit immediately with only start set as reachable set
% same goes for tStart = tFinal, which may occur in hybrid systems
if withinTol(params.tStart,params.tFinal) || ...
        ( isfield(options,'specification') && ~isempty(options.specification) ...
        && strcmp(options.specification(1).type,'invariant') ...
        && representsa_(options.specification(1).set,'emptySet',eps))
    timePoint.set{1} = params.R0; timePoint.time{1} = params.tStart;
    timeInt = [];
    res = false;
else

    % initialize taylorLinSys helper property for algorithms below
    if isempty(linsys.taylor)
        linsys.taylor = taylorLinSys(linsys.A);
    end

    % decide which reach function to execute by options.linAlg
    if strcmp(options.linAlg,'adaptive')
        [timeInt,timePoint,res] = reach_adaptive(linsys,params,options);
    else
        % all below, const. time step sizes
        if strcmp(options.linAlg,'standard')
            [timeInt,timePoint,res] = reach_standard(linsys,params,options);
        elseif strcmp(options.linAlg,'wrapping-free')
            [timeInt,timePoint,res] = reach_wrappingfree(linsys,params,options);
        elseif strcmp(options.linAlg,'fromStart')
            [timeInt,timePoint,res] = reach_fromStart(linsys,params,options);
        elseif strcmp(options.linAlg,'decomp')
            [timeInt,timePoint,res] = reach_decomp(linsys,params,options);
        elseif strcmp(options.linAlg,'krylov')
            [timeInt,timePoint,res] = reach_krylov(linsys,params,options);
        end
        % error vector (initial set: no error; error not computed -> NaN)
        timePoint.error = [0; NaN(length(timePoint.set)-1,1)];
        if isfield(timeInt,'set')
            timeInt.error = NaN(length(timeInt.set),1);
        end
    end
end

% delete all helper variables
linsys.taylor = [];

% create object of class reachSet
R = reachSet.initReachSet(timePoint,timeInt);

% check temporal logic specifications
if res && ~isempty(specLogic)
    res = check(specLogic,R);
    options.specification = spec;
end

% ------------------------------ END OF CODE ------------------------------
