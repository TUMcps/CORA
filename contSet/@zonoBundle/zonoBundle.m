classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonoBundle  < contSet
% zonoBundle - object constructor for zonotope bundles
%
% Description:
%    This class represents zonotope bundle defined as
%    \bigcap_j=1^k {c_j + \sum_{i=1}^p_j beta_i * g_j^(i) | beta_i \in [-1,1]},
%    i.e., the intersection of k zonotopes
%
% Syntax:  
%    obj = zonoBundle(list)
%
% Inputs:
%    list - cell-array list = {Z1,Z2,...} storing the zonotopes that
%           define the zonotope bundle
%
% Outputs:
%    obj - zonoBundle object
%
% Example:
%    Z1 = zonotope([1 3 0; 1 0 2]);
%    Z2 = zonotope([0 2 2; 0 2 -2]);
%    zB = zonoBundle({Z1,Z2});
%
%    figure; hold on;
%    plot(zB,[1,2],'FaceColor','r');
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:       Matthias Althoff
% Written:      09-November-2010
% Last update:  14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    Z;
    parallelSets;
end
    
methods
    %class constructor
    function obj = zonoBundle(varargin)
        
        if nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs', 1));
        end

        % parse input
        Z = setDefaultValues({{}}, varargin);

        if isa(Z, 'zonoBundle')
            % copy contructor
            obj = Z;
            return
        end

        inputArgsCheck({{Z, 'att', 'cell'}})

        % check dimensions
        if ~all(cellfun(@(x) isa(x,'zonotope'), Z,'UniformOutput',true))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'First input argument has to be a list of zonotope objects.'));
        end
        % all zonotopes have to be of the same dimension
        if any(diff(cellfun(@(x) dim(x), Z,'UniformOutput',true)))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Zonotopes have to be embedded in the same affine space.'));
        end

        % assign properties
        obj.Z = Z;
        obj.parallelSets = length(Z);
        
        % set parent object properties
        if isempty(obj.Z)
            obj.dimension = 0;
        else
            obj.dimension = size(obj.Z{1}.Z,1);
        end
    end
         
    % methods in seperate files
    zB = and(zB,S) % intersection
    zB = cartProd(zB,S) % Cartesian product
    c = center(zB) % center (only approximation)
    cPZ = conPolyZono(zB) % conversion to conPolyZono object
    res = contains(zB,S,varargin) % containment check
    zB = convHull(zB,varargin) % convex hull
    cZ = conZonotope(zB) % conversion to conZonotope object
    n = dim(zB) % dimension
    zB = enclose(zB,varargin) % enclose zonotope bundle and affine transformation
    zB = encloseTight(zB1,zB2,W) % enclose zonotope bundle and affine transformation
    zB = enlarge(zB,factor) % enlarge extensions
    I = interval(zB) % conversion to interval object
    res = isempty(zB) % emptyness check
    res = isequal(zB1,zB2) % equality check
    res = isFullDim(zB) % full-dimensionality check
    res = isIntersecting(zB,S,varargin) % intersection check
    P = mptPolytope(zB) % conversion to polytope
    zB = mtimes(factor1,factor2) % overloaded * operator
    zB = or(zB1, varargin) % union
    zB = plus(summand1,summand2) % overloaded + operator
    pZ = polyZonotope(zB) % conversion to polyZonotope object
    zB = project(zB,dims) % projection onto subspace
    zB = quadMap(zB,Q) % quadratic map
    p = randPoint(zB,varargin) % random point
    zB = reduce(zB,option,varargin) % order reduction
    zB = reduceCombined(zB,option,varargin) % order reduction
    zB = replace(zB,index,Z) % replace individual zonotope in bundle
    zB = shrink(zB,filterLength) % shrink extensions
    Zsplit = split(zB,varargin) % split
    [val,x] = supportFunc(zB,dir,varargin) % evaluate support function
    V = vertices(zB) % compute vertices
    res = volume(zB) % compute volume
    Z = zonotope(zB) % conversion to zonotope object
        
    %display functions
    han = plot(zB,varargin) % plot set
    display(zB) % display to console

end

methods (Static = true)
    zB = generateRandom(varargin) % generate random zonotope bundle
end


end

%------------- END OF CODE --------------