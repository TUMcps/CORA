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

% Authors:       Matthias Althoff
% Written:       09-November-2010
% Last update:   14-December-2022 (TL, property check in inputArgsCheck)
%                29-March-2023 (TL, optimized constructor)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------


properties (SetAccess = private, GetAccess = public)
    Z;              % list of zonotopes

    % internally-set properties
    parallelSets;   % number of zonotopes
end
    
methods
    % class constructor
    function obj = zonoBundle(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'zonoBundle')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        Z = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(Z,nargin);

        % 4. compute internal properties
        parallelSets = length(Z);

        % 5. assign properties
        obj.Z = Z;
        obj.parallelSets = parallelSets;

    end
         
    % methods in seperate files
    c = center(zB) % center (only approximation)
    cPZ = conPolyZono(zB) % conversion to conPolyZono object
    zB = convHull(zB,varargin) % convex hull
    cZ = conZonotope(zB) % conversion to conZonotope object
    n = dim(zB) % dimension
    zB = enclose(zB,varargin) % enclose zonotope bundle and affine transformation
    zB = encloseTight(zB1,zB2,W) % enclose zonotope bundle and affine transformation
    zB = enlarge(zB,factor) % enlarge extensions
    I = interval(zB) % conversion to interval object
    res = isequal(zB1,zB2) % equality check
    res = isFullDim(zB) % full-dimensionality check
    P = polytope(zB) % conversion to polytope
    zB = mtimes(factor1,factor2) % overloaded * operator
    zB = or(zB1, varargin) % union
    zB = plus(summand1,summand2) % overloaded + operator
    pZ = polyZonotope(zB) % conversion to polyZonotope object
    zB = project(zB,dims) % projection onto subspace
    zB = quadMap(zB,Q) % quadratic map
    zB = reduce(zB,option,varargin) % order reduction
    zB = reduceCombined(zB,option,varargin) % order reduction
    zB = replace(zB,index,Z) % replace individual zonotope in bundle
    res = representsa_(zB,type,tol,varargin) % comparsion to other set representations
    zB = shrink(zB,filterLength) % shrink extensions
    Zsplit = split(zB,varargin) % split
    Z = zonotope(zB) % conversion to zonotope object
        
    %display functions
    han = plot(zB,varargin) % plot set
    display(zB) % display to console

end

methods (Static = true)
    zB = generateRandom(varargin) % generate random zonotope bundle
    zB = empty(n) % instantiates an empty zonotope bundle
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 1
        throw(CORAerror('CORA:tooManyInputArgs',1));
    end

    % set default values
    Z = setDefaultValues({{}},varargin);

end

function aux_checkInputArgs(Z,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({{Z, 'att', 'cell'}})
    
        % check if zonotopes
        if ~all(cellfun(@(x) isa(x,'zonotope'), Z,'UniformOutput',true))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'First input argument has to be a list of zonotope objects.'));
        end
        % all zonotopes have to be of the same dimension
        if any(diff(cellfun(@(x) dim(x), Z,'UniformOutput',true)))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Zonotopes have to be embedded in the same affine space.'));
        end
        
    end

end

% ------------------------------ END OF CODE ------------------------------
