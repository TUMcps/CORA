function [res,cert,scaling] = contains_(Z,S,method,tol,maxEval,certToggle,scalingToggle)
% contains_ - determines if a zonotope contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(Z,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object or single point or matrix of points
%    method - method used for the containment check.
%       The available options are:
%           - 'exact': Checks for containment by using either
%               'exact:venum' or 'exact:polymax', depending on the number
%               of generators of Z and the object S.
%           - 'approx': Checks for containment using 'approx:st' (see  
%               below) if S is a zonotope, or any approximative method 
%               available otherwise.
%           - 'exact:venum': Checks for containment by enumerating all 
%               vertices of S (see Algorithm 1 in [2]).
%           - 'exact:polymax': Checks for containment by maximizing the
%               polyhedral norm w.r.t. Z over S (see Algorithm 2 in [2]).
%           - 'approx:st': Solves the containment problem using the
%               approximative method from [1]. If a solution using
%               'approx:st' returns that Z1 is contained in Z2, then this
%               is guaranteed to be the case. The runtime is polynomial
%               w.r.t. all inputs.
%           - 'approx:stDual': Solves the containment problem using the 
%               dual approximative method from [3]. Returns the same values
%               for res and scaling as 'approx:st', but cert can be more
%               precise.
%          For the next methods, note that if both certToggle and
%          scalingToggle are set to 'false', then res will be set to
%          'false' automatically, and the algorithms will not be executed.
%          This is because stochastic/optimization-based algorithms can not
%          confirm containment, so res = true can never happen. However, if
%          maxEval is set high enough, and res = false but cert = false,
%          one might conclude that with good probability, containment
%          holds.
%           - 'opt': Solves the containment problem via optimization
%               (see [2]) using the subroutine ga. If a solution
%               using 'opt' returns that Z1 is not contained in Z2, then
%               this is guaranteed to be the case. The runtime is
%               polynomial w.r.t. maxEval and the other inputs.
%           - 'sampling:primal': Solves the containment stochastically,
%               using the Shenmaier vertex sampling from [4].
%           - 'sampling:dual': Solves the containment stochastically, using
%               the Shenmaier halfspace sampling from [4].
%       The methods 'exact:venum' and 'exact:polymax' are only available if
%       S is a zonotope or a point/point cloud, and 'opt', 'approx:st', and
%       'approx:stDual' are only available if S is a zonotope.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - only, if 'opt', 'sampling:primal', or 'sampling:dual' is
%       used: Number of maximal function evaluations.
%    certToggle - if set to 'true', cert will be computed (see below).
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in Z, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in Z).
%           If res=true, then cert=true.
%           Note that computing this certification may marginally increase
%           the runtime.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z - center(Z)) + center(Z) contains S.
%           For the methods 'approx' and 'approx:st' this is an upper
%           bound, for 'opt', 'sampling:primal' and 'sampling:dual', this
%           number is a lower bound.
%           Note that computing this scaling factor may significantly 
%           increase the runtime.
%
% Note: For low dimensions or number of generators, and if S is a point
% cloud with a very large number of points, it may be beneficial to convert
% the zonotope to a polytope and call its containment operation
%
% Example: 
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
% 
%    contains(Z1,Z2)
%    contains(Z1,Z3)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%    
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%    [2] A. Kulmburg, M. Althoff.: On the co-NP-Completeness of the
%        Zonotope Containment Problem, European Journal of Control 2021
%    [3] A. Kulmburg, M. Althoff.: Hardness and Approximability of the
%        Containment Problem for Zonotopes and Ellipsotopes
%        (to appear)
%    [4] Kulmburg A., Brkan I., Althoff M.,: Search-based and Stochastic
%        Solutions to the Zonotope and Ellipsotope Containment Problems,
%        ECC 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains_, conZonotope/contains_

% Authors:       Matthias Althoff, Niklas Kochdumper, Adrian Kulmburg, Ivan Brkan
% Written:       07-May-2007 
% Last update:   06-April-2017
%                14-September-2019
%                19-November-2019 (NK, changed to header format)
%                01-July-2021 (AK, modified input parsing, implemented methods from [2])
%                22-July-2022 (MA, method st no longer requires YALMIP)
%                25-November-2022 (LS, method st using sparse matrices)
%                25-November-2022 (MW, rename 'contains')
%                05-February-2024 (AK, moved subfunctions, added sampling methods and stDual)
%                06-March-2024 (TL, check emptiness of zonotopes)
%                27-September-2024 (MW, remove halfspace call)
%                02-October-2024 (MW, point-in-zono, type decides LP/polytope)
%                15-January-2025 (TL, point-in-zono, tol is used for degeneracy check)
%                28-March-2025 (TL, buffer degenerate sets)
%                28-May-2025 (TL, quick check for representsa interval)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% No matter what, the algorithms will run better if we simplify Z as much
% as possible
Z = compact(Z);

% Quick checks for trivial cases ------------------------------------------

% represents a point?
[Z_isPoint, p] = representsa(Z, 'point');
if Z_isPoint
    if isnumeric(S)
        % point containment
        res = all(withinTol(S,p,tol));
        cert = true;
        if res
            scaling = 0;
        else
            scaling = Inf;
        end
        return
    else
        % check if S is also a point
        [S_isPoint, q] = representsa(S, 'point');
        if S_isPoint
            res = all(withinTol(p,q,tol));
            cert = true;
            if res
                scaling = 0;
            else
                scaling = Inf;
            end
        end
        return
    end
end

% represents an interval?
[Z_isInterval,I] = representsa_(Z,'interval',tol);
if Z_isInterval && any(startsWith(method,{'exact','approx'}))
    try
        [res,cert,scaling] = contains_(I,S,method,tol,maxEval,certToggle,scalingToggle);
    catch ME
        % check if a specific method was used
        if contains(method,':')
            % try with base method
            method = split(method,':');
            method = method{1};
            try
                 % retry with base method
                [res,cert,scaling] = contains_(I,S,method,tol,maxEval,certToggle,scalingToggle);
            catch ME2
                % not successfull, rethrow original exception
                rethrow(ME);
            end
        else
            % unable to fix automatically, rethrow exception
            rethrow(ME);
        end
    end
    return
end

% full check --------------------------------------------------------------

% check if full-dimensional
if ~isFullDim(Z,tol)
    % buffer degenerate sets slightly
    I = tol*interval(-ones(dim(Z),1),ones(dim(Z),1));
    Z = Z + I;
end

% point or point cloud in zonotope containment
if isnumeric(S)
    % Making sure the right algorithm is used
    if ~(strcmp(method, 'exact') || strcmp(method, 'exact:venum') || ...
            strcmp(method, 'exact:polymax') || strcmp(method, 'approx'))
        throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
    end

    [res, cert, scaling] = priv_zonotopeContainment_pointContainment(Z, S, method, tol, scalingToggle);
    return
end

% So, S is a contSet. Again, the simpler S is, the faster the algorithms
% are going to run:
S = compact(S);

% First, deal with trivial cases
if ~isBounded(S)
    % Unbounded -> not contained, since Z is always
    % bounded
    res = false;
    cert = true;
    scaling = Inf;
    return
elseif representsa(S, 'emptySet')
    % Empty -> always contained
    res = true;
    cert = true;
    scaling = 0;
    return
else
    try
        [isPoint,p] = representsa(S, 'point');
        if isPoint
            [res, cert, scaling] = priv_zonotopeContainment_pointContainment(Z,p,method,tol,scalingToggle);
            return
        end
    catch ME
        if strcmp('CORA:notSupported', ME.identifier) || strcmp('MATLAB:maxlhs',ME.identifier)
            % If the code above returns an error either because there are
            % too many outputs, or the operation is not supported, we
            % conclude that it is not implemented, and we do nothing
        else
            % In any other case, something went wrong. Relay that
            % information.
            rethrow(ME);
        end
    end
end


% Now that the trivial cases have been dealt with, we need to choose
% which algorithm to use.

switch method
    case {'exact','exact:venum','exact:polymax'} % Exact algorithms
        [res, cert, scaling] = aux_exactParser(Z, S, method, tol, maxEval, certToggle, scalingToggle);
    case {'approx','approx:st','approx:stDual'} % Approximative algorithms
        [res, cert, scaling] = aux_approxParser(Z, S, method, tol, maxEval, certToggle, scalingToggle);
    case {'sampling','sampling:primal','sampling:dual'} % Sampling algorithms
        [res, cert, scaling] = aux_samplingParser(Z, S, method, tol, maxEval, certToggle, scalingToggle);
    case 'opt' % Optimizer algorithms
        % Currently, for opt we only allow S to be a zonotope
        if ~isa(S, 'zonotope')
            throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
        end

        % Check if genetic algorithm optimization is available
        GOT_installed = true;
        try
            optimoptions('ga');
        catch
            GOT_installed = false;
        end

        if ~GOT_installed
            % show warning
            CORAwarning('CORA:solver', ...
                ['You have not installed the Global ' ...
                'Optimization Toolbox from MATLAB, and can '...
                'therefore not use the genetic algorithm for solving '...
                'the zonotope containment problem. '...
                'Alternatively, the DIRECT algorithm will '...
                'be used for now, but for improved results, '...
                'please install the Global Optimization Toolbox.']);
            [res,cert,scaling] = priv_zonotopeContainment_DIRECTMaximization(S, Z, tol, maxEval, scalingToggle);
        else
            % compute genetic
            [res,cert,scaling] = priv_zonotopeContainment_geneticMaximization(S, Z, tol, maxEval, scalingToggle);
        end
    otherwise
        throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
end

end

            
% Auxiliary functions -----------------------------------------------------

function [res, cert, scaling] = aux_exactParser(Z, S, method, tol, maxEval, certToggle, scalingToggle)
    % Chooses what algorithm to call upon to check for exact containment
    switch class(S)
        % If S is convex, there are two main options: Either S is a
        % polyhedron, or it is not.
        % In the first case, it makes sense to compute the vertices
        % of S, in the latter we have to compute the halfspaces of
        % Z.
        % Also, if the set is a zonotope, it will be practical to use the
        % zonotope-zonotope algorithms.
        case {'interval', 'zonotope'}
            % As of now, to respect backwards compatibility, the default
            % method for 'exact' is 'exact:venum' if S is an interval or
            % zonotope, 'exact:polymax' otherwise, except for dimensions
            % belowe 4
            S = zonotope(S);
            if strcmp(method, 'exact') && dim(S) >= 4
                method = 'exact:venum';
            else
                method = 'exact:polymax';
            end
        case {'conHyperplane', 'emptySet', 'fullspace', ...
                'halfspace', 'polytope', 'conZonotope', 'zonoBundle'}
            % Choose polymax, unless the user decides otherwise
            if strcmp(method, 'exact')
                method = 'exact:polymax';
            end

        case {'capsule', 'ellipsoid'}
            % Choose polymax, unless the user decides otherwise
            if strcmp(method, 'exact')
                method = 'exact:polymax';
            elseif strcmp(method, 'exact:venum')
                % Actually nevermind, venum is not an option here
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
            end

        otherwise
            throw(CORAerror('CORA:noExactAlg',Z,S));
    end

    % Now, it only remains to execute the corresponding algorithm
    switch method
        case 'exact:venum'
            if isa(S, 'zonotope')
                % Specialized function if S is a zonotope, as we can
                % exploit the fact that it can be written as a search-tree
                [res, cert, scaling] = priv_zonotopeContainment_vertexEnumeration(S, Z, tol, scalingToggle);
            else
                % All other cases are dealt with in the conZonotope
                % contains function
                [res, cert, scaling] = contains_(conZonotope(Z),S,'exact:venum',tol,maxEval,certToggle,scalingToggle);
            end
        case 'exact:polymax'
                % Send to polytops, the algorithm is implemented there
                % already
                [res, cert, scaling] = contains_(polytope(Z),S,'exact:polymax',tol,maxEval,certToggle,scalingToggle);
        otherwise
            throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
    end
end

function [res, cert, scaling] = aux_approxParser(Z, S, method, tol, maxEval, certToggle, scalingToggle)
    % Chooses what algorithm to call upon to check for approx containment

    switch class(S)
        % If S is convex, there are two main options: Either S is a
        % polyhedron, or it is not.
        % In the first case, we can use the methods from conZonotoep,
        % otherwise we first need to transform the set into a
        % polyhedron (right now, this is done by computing the interval
        % over-approximation)
        % Additionally, if S is a zonotope, we can use specific methods
        % from [1] and [3].
        case {'interval', 'zonotope'}
            S = zonotope(S);
            % If the user did not specifically choose st or stDual, the
            % main driver should be whether the user wants a certificate.
            % If yes, then stDual can actually be more accurate
            if certToggle && strcmp(method, 'approx')
                method = 'approx:stDual';
            elseif strcmp(method, 'approx')
                method = 'approx:st';
            end

        case {'conHyperplane', 'emptySet', 'fullspace', ...
                'halfspace', 'polytope', 'zonoBundle'}
            % Those are all the cases that can represent polyhedra
            % Filter out everything but 'approx'
            if ~strcmp(method, 'approx')
                % This is to avoid people using st here, as the algorithm
                % that is being used here is not quite st.
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
            end
        case {'capsule', 'ellipsoid'}
            % All other convex set representations
            if ~strcmp(method, 'approx')
                % This is to avoid people using st here, as the algorithm
                % that is being used here is not quite st.
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
            end
           
        otherwise
            % If S is not a convex set representation, our best chance is
            % to check for containment using the halfspace representation
            % of Z. But this is done in conZonotope/contains_, so we just
            % need to double-check that st and stDual are not used

           if ~strcmp(method, 'approx')
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
           end
    end

    % And now, it remains to execute the algorithm:
    switch method
        case 'approx'
            % This can only happen if S is not a zonotope
            % -> go to conZonotope/contains_
            [res, cert, scaling] = contains_(conZonotope(Z),S,method,tol,maxEval,certToggle,scalingToggle);
        case 'approx:st'
            [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrake(S, Z, tol, scalingToggle);
        case 'approx:stDual'
            [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrakeDual(S, Z, tol, scalingToggle);
        otherwise
            throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
    end
end

function [res, cert, scaling] = aux_samplingParser(Z, S, method, tol, maxEval, certToggle, scalingToggle)
    % Depending on the set type of S, we use a different algorithm
    switch class(S)
        case 'conZonotope'
            % The conZonotope case is a bit apart, since the algorithms do
            % not have the same precision bounds as in the other cases.
            % Consequently, we redirect to the conZonotope implementation,
            % so as to not confuse anybody.
            [res, cert, scaling] = contains_(conZonotope(Z), S,method,tol,maxEval,certToggle,scalingToggle);
        case 'zonotope'
            % The most straightforward case here. Use one of the dedicated
            % methods:
            if strcmp(method, 'sampling') || strcmp(method, 'sampling:primal')
                [res,cert,scaling] = priv_zonotopeContainment_zonoSampling(S, Z, tol, maxEval, scalingToggle);
            elseif strcmp(method, 'sampling:dual')
                [res,cert,scaling] = priv_zonotopeContainment_zonoSamplingDual(S, Z, tol, maxEval, scalingToggle);
            else
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
            end
        case 'ellipsoid'
            % This is an ellipsoid, so similar techniques as for zonotopes
            % can be used here. Still, since these methods are a bit
            % different (actually, even a bit simpler), we call upon
            % dedicated functions
            if strcmp(method, 'sampling') || strcmp(method, 'sampling:primal')
                [res,cert,scaling] = priv_zonotopeContainment_ellipsoidSampling(S, Z, tol, maxEval, scalingToggle);
            elseif strcmp(method, 'sampling:dual')
                [res,cert,scaling] = priv_zonotopeContainment_ellipsoidSamplingDual(S, Z, tol, maxEval, scalingToggle);
            else
                throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
            end
        otherwise
            throw(CORAerror('CORA:noSpecificAlg',method,Z,S));
    end

end


% ------------------------------ END OF CODE ------------------------------
