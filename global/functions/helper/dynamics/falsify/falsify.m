function [res,fals] = falsify(fun,params,varargin)
% falsify - find a falsifying trajectory for a black-box system
%
% Syntax:
%    [res,fals] = falsify(fun,params,spec)
%    [res,fals] = falsify(fun,params,options,spec)
%
% Inputs:
%    fun - function handle [t,x] = f(x0,u) for the black box simulation
%          function, where x0 is the initial state, u are the system 
%          inputs, and x and t are the states and time points of the
%          resulting trajectory
%    params - parameter defining the reachability problem
%       .R0: initial set (class contSet)
%       .U:  input set (class contSet)
%    options - options for falsification
%       .falsifyAlg: 'koopman'(default),'monteCarlo'
%       .nrConstInp: number of piecewise-constant input segments 
%                    (default: [], number is determined automatically)
%       .maxTime:    maximum computation time allocated for falsification 
%                    in seconds (default: 600)
%    spec - object of class specification (reach-avoid) or stl
%
% Outputs:
%    res - true/false whether falsification was successfull
%    fals - struct containing falsifying trajectory
%           .x0   ... point from initial set
%           .u    ... piecewise-constant input values
%           .tu   ... switching times of .u
%           .traj ... object of class trajectory
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       24-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check number of inputs
    narginchk(3,4);

    if nargin == 3
        spec = varargin{1}; options = [];
    else
        options = varargin{1}; spec = varargin{2};
    end
    
    % validate inputs
    if isempty(options)
        options = struct();
    end

    inputArgsCheck({ ...
        {fun, 'att', 'function_handle'}; ...
        {params, 'att', 'struct'}; ...
        {options, 'att', 'struct'}; ...
        {spec, 'att', {'specification','stl'}}; ...
    });

    if ~isa(spec,'specification')
        spec = specification(spec,'logic');
    end

    % set default algorithm parameters
    if ~isfield(params,'R0')
        params.R0 = interval(0);
    end

    if ~isfield(params,'U')
        params.U = interval(0);
    end

    % validate algorithm settings
    options = aux_checkOptions(options);

    % call the selected algorithm
    switch options.falsifyAlg
        case 'monteCarlo'
            [res,fals] = falsifyMonteCarlo(fun,params,options,spec);
        case 'koopman'
            [res,fals] = falsifyKoopman(fun,params,options,spec);
    end
end


% Auxiliary functions -----------------------------------------------------

function parsed = aux_checkOptions(options)
% check the algorithm settings provided by the user

    % default settings
    parsed.falsifyAlg = 'koopman';
    parsed.nrConstInp = 10;
    parsed.maxTime = 600;
    parsed.verbose = false;

    % check which options are provided by the user
    if ~isempty(options)
        
        if isfield(options,'falsifyAlg')
            parsed.falsifyAlg = options.falsifyAlg;
        end

        if isfield(options,'nrConstInp')
            parsed.nrConstInp = options.nrConstInp;
        end

        if isfield(options,'maxTime')
            parsed.maxTime = options.maxTime;
        end

        if isfield(options,'verbose')
            parsed.verbose = options.verbose;
        end
    end

    % check user defined settings
    if ~ismember(parsed.falsifyAlg,{'koopman','monteCarlo'})
        throw(CORAerror('CORA:wrongFieldValue','options.falsifyAlg','''koopman'' or ''monteCarlo'''));
    end

    if parsed.nrConstInp < 1 || mod(parsed.nrConstInp,1) ~= 0
        throw(CORAerror('CORA:wrongFieldValue','options.nrConstInp','integer > 0'));
    end

    if parsed.maxTime <= 0
        throw(CORAerror('CORA:wrongFieldValue','options.maxTime','double > 0'));
    end

    if ~islogical(parsed.verbose)
        throw(CORAerror('CORA:wrongFieldValue','options.verbose','boolean'));
    end

    % check if there are any redundant options specified
    redundantOptions(parsed,{'falsifyAlg','nrConstInp','maxTime','verbose'});
end

% ------------------------------ END OF CODE ------------------------------
