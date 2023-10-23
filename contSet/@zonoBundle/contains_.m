function res = contains_(zB,S,type,tol,varargin)
% contains_ - determines if a zonotope bundle contains a set or a point
%
% Syntax:
%    res = contains_(zB,S)
%    res = contains_(zB,S,type)
%    res = contains_(zB,S,type,tol)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object or single point
%    type - type of containment check:
%           'exact': any exact evaluation (see below, default)
%           'exact:zonotope': check containment for each zonotope
%           'exact:polytope: convert zB to polytope
%           'approx'
%
% Outputs:
%    res - true/false
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

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 16-March-2023 (MW, restructure, add more types)
%                27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% point or point cloud in zonotope bundle containment
if isnumeric(S)
    
    if strcmp(type,'exact:polytope')
        % conversion to polytope
        res = contains_(polytope(zB),S,'exact',tol);

    elseif strcmp(type,'exact:zonotope')
        % check every zonotope individually
        res = true(1,size(S,2));
        for i=1:zB.parallelSets
            res = res & contains_(zB.Z{i},S,'exact',tol);
        end

    else % all other types, including approx
        % conversion to constrained zonotope
        res = contains_(conZonotope(zB),S,'exact',tol);

    end
    
% capsule/ellipsoid in zonotope bundle containment
elseif isa(S,'capsule') || isa(S,'ellipsoid')
    % same algorithm for all types
    
    if strcmp(type,'exact:zonotope')
        % check every zonotope individually
        res = true(1,size(S,2));
        for i=1:zB.parallelSets
            res = res & contains_(zB.Z{i},S,'exact',tol);
        end

    else
        % default: conversion to polytope
        P = polytope(zB);
        res = contains_(P,S,'exact',tol);

    end

elseif isa(S,'interval')

    if strcmp(type,'exact:polytope')
        res = contains_(polytope(zB),S,'exact',tol);

    elseif strcmp(type,'exact:zonotope')
        % check every zonotope individually
        res = true(1,size(S,2));
        for i=1:zB.parallelSets
            res = res & contains_(zB.Z{i},S,'exact',tol);
        end
    
    else % all other types: check every vertex
        res = all(contains_(zB,vertices(S),'exact',tol));

    end

elseif isa(S,'taylm') || isa(S,'polyZonotope')

    % no exact algorithm
    if contains(type,'exact')
        throw(CORAerror('CORA:noExactAlg',zB,S));
    end

    % approx method via conversion to polytope
    res = contains_(polytope(zB),S,'exact',tol);

else % all other classes

    if contains(type,'exact')
        res = contains_(polytope(zB),S,type,tol);
    elseif strcmp(type,'approx')
        res = contains_(conZonotope(zB),S,type,tol);
    end

end

% ------------------------------ END OF CODE ------------------------------
