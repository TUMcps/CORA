function [res,cert,scaling] = contains_(E,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if an ellipsoid contains a set or a point
%
% Syntax:
%    res = contains_(E,S)
%    res = contains_(E,S,mode)
%
% Inputs:
%    E - ellipsoid object 
%    S - contSet object or single point
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of C
%       will be detected as lying in C, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in E, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in E).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(E - E.center) + E.center contains S.
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2));
%    Z = zonotope([0 1 0;0 1 1]);
%
%    contains(E1,E2)
%    contains(E1,Z)
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%    plot(Z,[1,2],'r');
%
% References:
%    [1] Yildirim, E.A., 2006. On the minimum volume covering ellipsoid of
%        of ellipsoids. SIAM Journal on Optimization, 17(3), pp.621-641.     
%    [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%    [3] Kulmburg, A, Schäfer, L, Althoff, A, 2024. Approximability of the
%        Containment Problem for Zonotopes and Ellipsotopes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Victor Gassmann, Niklas Kochdumper, Adrian Kulmburg, Mark Wetzlinger
% Written:       15-October-2019 
% Last update:   21-November-2019 (NK, extend to other sets)
%                09-March-2021 (included tolerance for q comparison)
%                17-March-2021 (error handling)
%                06-July-2021 (AK, merged containsPoint to in)
%                04-July-2022 (VG, handle class array cases)
%                25-November-2022 (MW, rename 'contains')
%                25-April-2023 (VG, add method for capsule)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

cert = 1;
scaling = 0;

% Check out trivial cases
[ell_isPoint, p] = representsa(E, 'point');
if ell_isPoint
    if isnumeric(S)
        res = max(max(abs(S-p))) <= tol;
        cert = true;
        if res
            scaling = 0;
        else
            scaling = Inf;
        end

    else % S is not numeric
        [S_isPoint, q] = representsa(S, 'point');
        if S_isPoint
            res = all(q==p);
            cert = true;
            if res
                scaling = 0;
            else
                scaling = Inf;
            end
        elseif representsa(S, 'emptySet')
            res = true;
            cert = true;
            scaling = 0;
        else % S is not numeric and not empty
            res = false;
            cert = true;
            scaling = inf;
        end
    end
    return
end
% E is not a point

if representsa(E, 'emptySet')
    if isnumeric(S)
        if isempty(S)
            res = true;
            scaling = 0;
            cert = true;
        else
            res = false;
            scaling = Inf;
            cert = true;
        end
    else
        if representsa(S, 'emptySet')
            res = true;
            scaling = 0;
            cert = true;
        else
            res = false;
            scaling = Inf;
            cert = true;
        end
    end
    return
end
% E is not empty

if isnumeric(S)
    [res,cert,scaling] = priv_containsPoint(E,S,tol);
    return
end

if ~isBounded(S)
    % Unbounded -> not contained, since E is always
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
            [res, cert, scaling] = aux_containsPoint(E,p,method,tol,scalingToggle);
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

% containment check


if isa(S,'ellipsoid')
    [res,cert,scaling] = priv_containsEllipsoid(E,S,tol);
    return
end

if isa(S,'capsule')
    % check if balls at both ends of capsule are contained
    E1 = ellipsoid(S.r^2*eye(dim(S)), S.c+S.g);
    E2 = ellipsoid(S.r^2*eye(dim(S)), S.c-S.g);

    [res1, cert1, scaling1] = priv_containsEllipsoid(E,E1,tol);
    [res2, cert2, scaling2] = priv_containsEllipsoid(E,E2,tol);
    res = res1 & res2;
    cert = true;
    scaling = max([scaling1 scaling2]);
    return
end

if strcmp(method, 'exact')
    if isa(S, 'conZonotope') || isa(S, 'interval') || ...
            isa(S, 'polytope') || isa(S, 'zonoBundle')
        % check if all vertices of the set are contained
        [res,cert,scaling] = priv_containsPoint(E,vertices(S),tol);
        res = all(res);
        cert = true;
        scaling = max(scaling);
    elseif isa(S, 'zonotope')
        % For zonotopes, we can leverage the symmetry for a better vertex
        % enumeration, at least for dimensions >3. For some reason,
        % computing the vertex representation in dimension 2 is still
        % faster:
        if dim(S) <= 2
            [res,cert,scaling] = priv_containsPoint(E,vertices(S),tol);
            res = all(res);
            cert = true;
            scaling = max(scaling);
        else
            [res,cert,scaling] = priv_venumZonotope(E,S,tol,scalingToggle);
        end
    else % throw error
        throw(CORAerror('CORA:noExactAlg',E,S));
    end
elseif strcmp(method, 'approx')
    % compute approx algorithms
    if isa(S,'zonotope') && strcmp(method,'approx')
        [res,cert,scaling] = aux_symmetricGrothendieck(E,S,tol,certToggle);
        return
    elseif ismethod(S,'zonotope')
        [res,cert,scaling] = contains_(E,zonotope(S),method,0,certToggle,scalingToggle);
        cert = res;
    else
        throw(CORAerror('CORA:noExactAlg',E,S));
    end
else
    throw(CORAerror('CORA:noSpecificAlg',E,S,method));
end
    
end


% Auxiliary functions -----------------------------------------------------

function [res, cert, scaling] = aux_symmetricGrothendieck(E,Z,tol,certToggle)

n = dim(Z);
m = size(Z.generators, 2);

% First, we deal with the centers; Z < E if and only if Z' < E-center(E),
% where Z' is the zonotope with generator matrix [G center(Z)-center(E)]:

Z = zonotope(zeros([n 1]), [Z.generators Z.center-center(E)]);
E = E - center(E);

G = Z.generators;

% We now deal with the case where E is degenerate:
if ~isFullDim(E)
    % Z < E can only happen if Z lies in the same subspace as E, so we need
    % to check that every generator of Z lies, up to some scaling, in E.
    % The fastest way to do this is via the ellipsoid norm: if it is equal
    % to Inf for some generator, this means that the generator in question
    % does not lie in the same subspace:
    for i=1:m
        if E.ellipsoidNorm(G(:,i)) == Inf
            res = false;
            cert = false;
            scaling = Inf;
            return
        end
    end
    
    % So, now we know that Z lies in the same subspace as E. Let us rotate
    % E using the svd in such a way, that the axes of E are aligned with
    % the canonical ONB:
    [T,~,~] = svd(E.Q);
    E = T'*E;
    % We need to rotate Z too:
    Z = T'*Z;
    
    % We can now remove the last coordinates
    r = rank(E.Q);
    E = ellipsoid(E.Q(1:r,1:r), zeros([r 1]));
    G = Z.generators;
    Z = zonotope(zeros([r 1]), G(1:r,:));
end

% We can now assume that E is non-degenerate, which means that Q is
% invertible.
% We now use the method described in [3]
G = Z.generators;
m = size(G,2);
X = sdpvar(m, m);
lambda = sdpvar(1,1);
constraints = [X >= 0];
for i=1:m
   constraints = [constraints X(i,i)==1]; 
end

Qinv = inv(E.Q);

cost = -trace(G'*Qinv*G*X);
    
    
options = sdpsettings('solver', 'sedumi', 'verbose',0, 'allownonconvex',0);

% solve PSD programming problem
yalmipOptimizer = optimizer(constraints, cost, options,[],{X});

optimizationResults = yalmipOptimizer();

sol = optimizationResults;

scaling = sqrt(abs(trace(G'*Qinv*G*sol)));
if scaling <= 1+tol
   res = true;
   cert = true;
else
   res = false;
   if scaling > pi/2
       cert = true;
   else
       cert = false;
   end
end

end


% ------------------------------ END OF CODE ------------------------------
