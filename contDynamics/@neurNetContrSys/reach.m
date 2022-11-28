function [R, res] = reach(obj, params, options, evParams, varargin)
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
%    options - options for the computation of reachable sets
%    evParams - parameters for neural network evaluation
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

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      17-September-2021
% Last update:  05-April-2022 (TL)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

options = aux_parseSettings(obj, params, options);

% adapt specification to extended system dimension
spec = [];
if nargin >= 5
    spec = aux_adaptSpecification(varargin{1}, obj.dim, obj.nn.neurons_out);
end

% loop over the whole time horizon
tVec = options.tStart:obj.dt:options.tFinal;
if tVec(end) ~= options.tFinal
    tVec = [tVec, options.tFinal];
end
X = options.R0;
R = [];

for i = 1:length(tVec) - 1

    % get control input = ouput of neural network
    U = obj.nn.evaluate(X, evParams);

    % update initial set and time
    if isa(X, 'zonotope')
        Z = [X.Z, zeros(size(X.Z, 1), size(U.Z, 2)-size(X.Z, 2)); U.Z];
        options.R0 = zonotope(Z);
    elseif isa(X, 'polyZonotope')
        if ~isempty(X.Grest)
            diffXU = max(size(X.Grest, 2)-size(U.Grest, 2), 0);
            diffUX = max(size(U.Grest, 2)-size(X.Grest, 2), 0);
            Grest = [X.Grest, zeros(size(X.Grest, 1), diffUX); ...
                zeros(size(U.Grest, 1), diffXU), U.Grest];
        else
            Grest = [zeros(dim(X), size(U.Grest, 2)); U.Grest];
        end

        if all(size(X.expMat) == size(U.expMat)) ...
                && all(all(X.expMat == U.expMat)) ... 
                && all(size(X.id) == size(U.id)) ...
                && all(X.id == U.id)
            options.R0 = polyZonotope([X.c; U.c], [X.G; U.G], Grest, ...
                X.expMat, X.id);
        else
            ids = unique([X.id; U.id]);
            
            c = [X.c; U.c];
            G = blkdiag(X.G, U.G);

            expMat = zeros(length(ids), size(G, 2));
            expMat(ismember(ids, X.id), 1:size(X.expMat, 2)) = X.expMat;
            expMat(ismember(ids, U.id), end-size(U.expMat, 2)+1:end) = U.expMat;

            options.R0 = polyZonotope(c, G, Grest, expMat, ids);
        end
    end
    options.tStart = tVec(i);
    options.tFinal = tVec(i+1);

    % compute reachable set for the controlled system
    [Rtemp, res] = aux_reachability(obj.sys, options, spec);

    % store reachable set
    R = add(R, project(Rtemp, 1:obj.dim));
    X = project(Rtemp.timePoint.set{end}, 1:obj.nn.neurons_in);
    X = X.replaceId(1:length(X.id));

    % terminate if specifications are violated
    if ~res
        return;
    end
end
end


% Auxiliary Functions -----------------------------------------------------

function spec = aux_adaptSpecification(spec, n, m)
% adapt the specifications to the extended system dimensions that include
% the control input as additional states

for i = 1:length(spec)
    poly = mptPolytope(spec(i).set);
    poly = projectHighDim(poly, n+m, 1:n);
    spec(i) = specification(poly, spec(i).type, spec(i).time);
end
end

function options = aux_parseSettings(obj, params, options)
% initilize the algorithm settings

% check if the algorithm settings provided by the user are correct
params.R0 = cartProd(params.R0, zeros(obj.nn.neurons_out, 1));
options = validateOptions(obj.sys, mfilename, params, options);
options.R0 = project(options.R0, 1:obj.nn.neurons_in);

% obtain factors for initial state and input solution time step
r = options.timeStep;
for i = 1:(options.taylorTerms + 1)
    options.factor(i) = r^(i) / factorial(i);
end

% initialize time-varying inputs
if ~isfield(options, 'uTransVec')
    temp = options.tStart:options.timeStep:options.tFinal;
    options.uTransVec = repmat(options.uTrans, [1, length(temp)]);
end

% check if splitting is turned off
if ~all(isinf(options.maxError))
    throw(CORAerror('CORA:notSupported',...
        ['Splitting reachable sets is not supported for neural', ...
        ' network controlled systems!']));
end

% pre-compute derivatives
derivatives(obj.sys, options);
end

function [R, res] = aux_reachability(obj, options, spec)
% compute the reachable set of the controlled system (without performing
% option checks)

res = true;
tVec = options.tStart:options.timeStep:options.tFinal;

% initialize cell-arrays that store the reachable set
Rint.set = cell(length(tVec)-1, 1);
Rint.time = Rint.set;
Rpoint.set = cell(length(tVec), 1);
Rpoint.time = Rpoint.set;
Rpoint.set{1} = options.R0;
Rpoint.time{1} = options.tStart;

% loop over all reachability steps
for i = 1:length(tVec) - 1

    % compute reachable set for one reachability step
    try
        options.uTrans = options.uTransVec(:, i);
        if i == 1
            [Rnext, options] = initReach(obj, options.R0, options);
        else
            [Rnext, options] = post(obj, Rnext, options);
        end
    catch ME
        R = aux_constructReachSet(Rpoint, Rint, i-1);

        if true || strcmp(ME.identifier, 'reach:setexplosion') || strcmp(ME.identifier, 'CORA:reachSetExplosion')
            % display information to user
            fprintf("\n");
            disp(ME.message);
            disp("  Step "+i+" at time t="+options.tStart);
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

%------------- END OF CODE --------------