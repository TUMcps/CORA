function [res,cert,scaling] = contains_(zB,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a zonotope bundle contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(zB,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object or single point
%    method - method used for the containment check.
%           'exact': any exact evaluation (see below, default)
%           'exact:zonotope': check containment for each zonotope
%           'exact:polytope': convert zB to polytope
%           'approx'
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       zB will be detected as lying in zB, which can be useful to
%       counteract errors originating from floating point errors.
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
%           guaranteed to not be contained in zB, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in zB).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(zB - zB.center) + zB.center contains S.
%
% Example: 
%    I1 = interval([0;-1],[2;2]);
%    I2 = I1 + [2;0];
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    contains(zB,I1)
%    contains(zB,I2)
%
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(I1,[1,2],'g');
%    
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(I2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains_, conZonotope/contains_, zonotope/contains_

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Adrian Kulmburg
% Written:       19-November-2019
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
%                16-January-2025 (AK, adding scaling and cert)
% Last revision: 16-March-2023 (MW, restructure, add more types)
%                27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% Deal with trivial cases first
if representsa(S, 'emptySet', tol)
    % Empty set is always contained
    res = true;
    cert = true;
    scaling = 0;
    return
elseif representsa(S, 'fullspace', tol)
    % Fullspace is never contained, since a zonotope bundle is compact
    res = false;
    cert = true;
    scaling = inf;
    return
end

% check emptySet
if representsa(zB, 'emptySet')
    if isnumeric(S)
        % If S is numeric, check if S is empty manually
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
        % If S is a contSet object, check whether it is empty; if not, it
        % can not possibly be contained
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

% point or point cloud in zonotope bundle containment
if isnumeric(S)
    if strcmp(method,'exact')
        [res,cert,scaling] = contains_(conZonotope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);
    elseif strcmp(method,'exact:polytope')
        % conversion to polytope
        [res,cert,scaling] = contains_(polytope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);

    elseif strcmp(method,'exact:zonotope')
        % check every zonotope individually
        res = true([1 size(S,2)]);
        cert = true([1 size(S,2)]);
        scaling = zeros([1 size(S,2)]);
        for i=1:zB.parallelSets
            [res_part, ~, scaling_part] = contains_(zB.Z{i},S,'exact',tol,maxEval,certToggle,scalingToggle);
            res = res & res_part;
            scaling = max([scaling;scaling_part]);
        end

    else % all other types, including approx
        % conversion to constrained zonotope
        [res,cert,scaling] = contains_(conZonotope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);
    end
    
% capsule/ellipsoid in zonotope bundle containment
elseif isa(S,'capsule') || isa(S,'ellipsoid')
    % same algorithm for all types
    
    if strcmp(method,'exact:zonotope')
        % check every zonotope individually
        res = true;
        cert = true;
        scaling = 0;
        for i=1:zB.parallelSets
            [res_part, ~, scaling_part] = contains_(zB.Z{i},S,'exact',tol,maxEval,certToggle,scalingToggle);
            res = res & res_part;
            scaling = max([scaling scaling_part]);
        end

    elseif contains(method, 'exact')
        % default: conversion to polytope
        P = polytope(zB);
        [res, cert, scaling] = contains_(P,S,'exact',tol,maxEval,certToggle,scalingToggle);
    else
        % approx method
        P = polytope(zB);
        [res, cert, scaling] = contains_(P,S,'approx',tol,maxEval,certToggle,scalingToggle);
    end

elseif isa(S,'interval')

    if strcmp(method,'exact:polytope')
        [res, cert, scaling] = contains_(polytope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);

    elseif strcmp(method,'exact:zonotope')
        % check every zonotope individually
        res = true;
        cert = true;
        scaling = 0;
        for i=1:zB.parallelSets
            [res_part, ~, scaling_part] = contains_(zB.Z{i},S,'exact',tol,maxEval,certToggle,scalingToggle);
            res = res & res_part;
            scaling = max([scaling scaling_part]);
        end
    
    else % all other types: check every vertex
        % (this includes 'exact' and 'approx', but 'approx' reduces to
        % 'exact')
        [res, cert, scaling] = contains_(zB,vertices(S),'exact',tol,maxEval,certToggle,scalingToggle);
        res = all(res);
        cert = true;
        scaling = max(scaling);
    end

elseif isa(S,'taylm') || isa(S,'polyZonotope') || isa(S,'conPolyZono')

    % no exact algorithm
    if contains(method,'exact')
        throw(CORAerror('CORA:noExactAlg',zB,S));
    end

    % approx method via conversion to polytope
    [res,cert,scaling] = contains_(polytope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);

else % all other classes

    if contains(method,'exact')
        [res,cert,scaling] = contains_(polytope(zB),S,'exact',tol,maxEval,certToggle,scalingToggle);
    elseif strcmp(method,'approx')
        [res,cert,scaling] = contains_(conZonotope(zB),S,'approx',tol,maxEval,certToggle,scalingToggle);
    end

end

% ------------------------------ END OF CODE ------------------------------
