function [R, res] = reach(obj, params, options, varargin)
% reach - computes the reachable continuous set for the entire time horizon
%    of a continuous system
%
% Syntax:
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - neurNetContrSys object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets and network evaluation
%    spec - object of class specification
%
% Outputs:
%    R - object of class reachSet storing the computed reachable set
%    res  - true if specifications are satisfied, otherwise false
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   05-April-2022 (TL)
%                28-March-2023 (TL, parse input, clean up)
%                20-July-2023 (TL, bugfix: mismatch time horizon & sampling)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% parse input
narginchk(3,4);
spec = setDefaultValues({[]}, varargin);
inputArgsCheck({ ...
    {obj, 'att', 'neurNetContrSys'}; ...
    {params, 'att', 'struct'}; ...
    {options, 'att', 'struct'}; ...
    {spec, 'att', {'specification', 'numeric'}, {{}, {'empty'}}}; ...
})

% validate user input
[params,options] = aux_parseSettings(obj, params, options);

% adapt specification
if isnumeric(spec)
    spec = [];
else
    spec = aux_adaptSpecification(spec, obj.nrOfStates, obj.nn.neurons_out);
end

% init
tVec = params.tStart:obj.dt:params.tFinal;
if tVec(end) ~= params.tFinal
    % add tFinal if sampling time and time horizon don't match
    % resulting in a partial time step at the end.
    tVec(end+1) = params.tFinal;
end
X = params.R0;
R = [];
idxLayer = 1:numel(obj.nn.layers);

for i = 1:length(tVec) - 1
    % compute next time step
    params.tStart = tVec(i);
    params.tFinal = tVec(i+1);

    % get control input = ouput of neural network
    U = obj.nn.evaluate_(X, options, idxLayer);

    % append control input to initial set
    params.R0 = aux_appendU2X(X, U);

    % compute reachable set for the controlled system of the next step
    % with constant input U
    [R_i, res] = aux_reachability(obj.sys, params, options, spec);

    % store reachable set
    R_i = project(R_i, 1:obj.nrOfStates);
    R = add(R, R_i);
    X = R_i.timePoint.set{end};

    % terminate if specifications are violated
    if ~res
        return;
    end
end
end


% Auxiliary functions -----------------------------------------------------

function spec = aux_adaptSpecification(spec, n, m)
% adapt the specifications to the extended system dimensions that include
% the control input as additional states

for i = 1:length(spec)
    poly = polytope(spec(i).set);
    poly = lift_(poly, n+m, 1:n);
    spec(i) = specification(poly, spec(i).type, spec(i).time);
end
end

function [params,options] = aux_parseSettings(obj, params, options)

% initilize the algorithm settings ----------------------------------------

% check if the algorithm settings provided by the user are correct
% doing some corrections to reuse config file for nonlinearSys/reach
params.R0 = cartProd(params.R0, zeros(obj.nn.neurons_out, 1));
optionsnn = options.nn; % is checked via nnHelper.validateNNoptions
options = rmfield(options,'nn');
[params,options] = validateOptions(obj.sys, params, options, 'FunctionName', 'reach');
params.R0 = project(params.R0, 1:obj.nn.neurons_in);
options.nn = optionsnn;

% obtain factors for initial state and input solution time step
r = options.timeStep;
for i = 1:(options.taylorTerms + 1)
    options.factor(i) = r^(i) / factorial(i);
end

% initialize time-varying inputs
if ~isfield(params, 'uTransVec')
    tVec = params.tStart:options.timeStep:params.tFinal;
    if tVec(end) ~= params.tFinal
        % add tFinal if sampling time and time horizon don't match
        % resulting in a partial time step at the end.
        tVec(end+1) = params.tFinal;
    end
    params.uTransVec = repmat(params.uTrans, [1, length(tVec)]);
end

% check if splitting is turned off
if ~all(isinf(options.maxError))
    throw(CORAerror('CORA:notSupported',...
        ['Splitting reachable sets is not supported for neural', ...
        ' network controlled systems!']));
end

% parse nn options --------------------------------------------------------

if ~isfield(options, 'add_approx_error_to_GI')
    % to reduce computational overhead
    options.nn.add_approx_error_to_GI = true;
end

options = nnHelper.validateNNoptions(options);

% pre-compute derivatives -------------------------------------------------

derivatives(obj.sys, options);

end

function R0 = aux_appendU2X(X, U)
    if isa(X, 'zonotope')
        % TODO: only works if no order reduction is applied within nn
        G = [X.G, zeros(size(X.G, 1), size(U.G, 2)-size(X.G, 2)); U.G];
        R0 = zonotope([X.c;U.c],G);
    elseif isa(X, 'polyZonotope')
        if ~isempty(X.GI)
            if size(X.GI, 2) > size(U.GI, 2)
                throw(CORAerror("CORA:wrongValue", ...
                    ['Not all generations in X.GI are still present' ...
                    ' in U.GI. Might be due to order reductions ' ...
                    'with the nn evaluation'] ...
                    ))
            end

            % TODO: only works if no order reduction is applied within nn
            diffUX = max(size(U.GI, 2)-size(X.GI, 2), 0);
            GI = [X.GI, zeros(size(X.GI, 1), diffUX); U.GI];
        else
            GI = [zeros(dim(X), size(U.GI, 2)); U.GI];
        end

        % TODO: exactPlus using right ids
        if all(size(X.E) == size(U.E)) ...
                && all(all(X.E == U.E)) ... 
                && all(size(X.id) == size(U.id)) ...
                && all(X.id == U.id)
            R0 = polyZonotope([X.c; U.c], [X.G; U.G], GI, X.E, X.id);
        else
            ids = unique([X.id; U.id]);
            
            c = [X.c; U.c];
            G = blkdiag(X.G, U.G);

            E = zeros(length(ids), size(G, 2));
            E(ismember(ids, X.id), 1:size(X.E, 2)) = X.E;
            E(ismember(ids, U.id), end-size(U.E, 2)+1:end) = U.E;

            R0 = polyZonotope(c, G, GI, E, ids);
        end
    end
end

function [R, res] = aux_reachability(sys, params, options, spec)
% compute the reachable set of the controlled system (without performing
% option checks)

res = true;

% init time steps
tVec = params.tStart:options.timeStep:params.tFinal;
if tVec(end) ~= params.tFinal
    % add tFinal if sampling time and intermediate time horizon don't match
    % resulting in a partial time step at the end.
    tVec(end+1) = params.tFinal;
end

% initialize cell-arrays that store the reachable set
Rint.set = cell(length(tVec)-1, 1);
Rint.time = Rint.set;
Rpoint.set = cell(length(tVec), 1);
Rpoint.time = Rpoint.set;
Rpoint.set{1} = params.R0;
Rpoint.time{1} = params.tStart;

% loop over all reachability steps
for i = 1:length(tVec) - 1
    options.timeStep = tVec(i+1) - tVec(i);

    % compute reachable set for one reachability step
    try
        params.uTrans = params.uTransVec(:, i);
        if i == 1
            [Rnext, options] = initReach(sys, params.R0, params, options);
            % initReach does not reduce
            Rnext.ti{1} = reduce(Rnext.ti{1},...
                options.reductionTechnique,options.zonotopeOrder);
            Rnext.tp{1}.set = reduce(Rnext.tp{1}.set,...
                options.reductionTechnique,options.zonotopeOrder);
        else
            [Rnext, options] = post(sys, Rnext, params, options);
        end
    catch ME
        R = aux_constructReachSet(Rpoint, Rint, i-1);

        if true || strcmp(ME.identifier, 'reach:setexplosion') ...
            || strcmp(ME.identifier, 'CORA:reachSetExplosion')
            % display information to user
            fprintf("\n");
            disp(ME.message);
            disp("  Step "+i+" at time t="+params.tStart);
            disp("The reachable sets until the current step are returned.");
            fprintf("\n");
            res = false;

        else
            % any other run-time error: report information
            rethrow(ME);
        end
        return;
    end

    % save reachable set
    Rint.set{i} = Rnext.ti{1};
    Rpoint.set{i+1} = Rnext.tp{1}.set;
    Rint.time{i} = interval(tVec(i), tVec(i+1));
    Rpoint.time{i+1} = tVec(i+1);

    % check specification
    if ~isempty(spec)
        res = check(spec, Rnext.ti, Rint.time{i});
        if ~res
            R = aux_constructReachSet(Rpoint, Rint, i);
            return;
        end
    end
end

% create resulting reachSet object
R = reachSet(Rpoint, Rint);

end

function R = aux_constructReachSet(Rpoint, Rint, i)
% construct a reachable set object

Rpoint.set = Rpoint.set(1:i+1);
Rint.set = Rint.set(1:i);
Rpoint.time = Rpoint.time(1:i+1);
Rint.time = Rint.time(1:i);

R = reachSet(Rpoint, Rint);

end

% ------------------------------ END OF CODE ------------------------------
