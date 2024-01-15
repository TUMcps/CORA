classdef (InferiorClasses = {?mp}) interval < contSet
% interval - object constructor for real-valued intervals 
%
% Description:
%    This class represents interval objects defined as
%    {x | a_i <= x <= b_i, \forall i = 1,...,n}.
%
% Syntax:
%    obj = interval(I)
%    obj = interval(a)
%    obj = interval(a,b)
%
% Inputs:
%    I - interval object
%    a - lower limit
%    b - upper limit
%
% Outputs:
%    obj - generated interval object
%
% Example:
%    a = [1;-1];
%    b = [2;3];
%    I = interval(a,b);
%    plot(I,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       19-June-2015
% Last update:   18-November-2015
%                26-January-2016
%                15-July-2017 (NK)
%                01-May-2020 (MW, delete redundant if-else)
%                20-March-2021 (MW, error messages)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                29-March-2023 (TL, optimized constructor)
%                08-December-2023 (MW, handle [-Inf,-Inf] / [Inf,Inf] case)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    inf;    % lower bound
    sup;    % upper bound
end

methods

    % class constructor
    function obj = interval(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'interval')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [lb,ub] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(lb,ub,nargin);

        % 4. compute properties (deal with corner cases)
        [lb,ub] = aux_computeProperties(lb,ub);

        % 5. assign properties
        obj.inf = lb;
        obj.sup = ub;

    end
    

    function ind = end(obj,k,n)
    % overloads the end operator for referencing elements, e.g. I(end,2),
        ind = size(obj,k);
    end
    
    % methods in seperate files
    res = abs(I) % absolute value function
    I = acos(I) % inverse cosine function
    I = acosh(I) % inverse hyperbolic cosine function
    I = asin(I) % inverse sine function
    I = asinh(I) % inverse hyperbolic sine function
    I = atan(I) % inverse tangent function
    I = atanh(I) % inverse hyperbolic tangent function
    C = capsule(I) % conversion to capsule object
    c = center(I) % center of interval
    cPZ = conPolyZono(I) % conversion to conPolyZono object
    res = convHull(I,varargin) % convex hull
    cZ = conZonotope(I) % conversion to conZonotope object
    res = cos(I) % cosine function
    I = cosh(I) % hyperbolic cosine function
    I = ctranspose(I) % overloaded ' operator
    res = diag(I,varargin) % overloaded diag-function
    n = dim(I) % dimension of interval
    E = ellipsoid(I) % conversion to ellipsoid object
    I = enlarge(I,factor) % enlargement by factor
    res = eq(I1,I2) % equality check
    I = exp(I) % overloaded exp-function
    p = gridPoints(I,segments) % generate grid points
    I = horzcat(varargin) % overloaded horizontal concatenation
    res = infimum(I) % read lower limit
    res = isequal(I1,I2,varargin) % equal objects check
    res = isFullDim(I) % full dimensionality check
    res = isscalar(I) % one-dimensionality check
    res = issparse(I) % issparse
    res = le(I1,I2) % subseteq check
    l = length(I) % largest dimension of interval
    I = lift_(I,N,proj) % lift to a high-dimensional space
    I = log(I) % logarithm function
    res = lt(I1,I2) % subset check
    I = max(I,Y,varargin) % maximum
    I = minkDiff(I,S,varargin) % Minkowski difference
    I = min(I,Y,varargin) % minimum 
    res = minus(minuend,subtrahend) % overloaded - operator (binary)
    res = mpower(base,exponent) % overloaded ^ operator
    P = polytope(I) % conversion to polytope object
    res = mrdivide(numerator,denominator) % overloaded / operator
    res = mtimes(factor1,factor2) % overloaded * operator
    res = ne(I1,I2) % overloaded ~= operator
    res = or(I,S) % union
    dzNew = partition(I, splits) % partition into subintervals
    han = plot(I,varargin) % plot
    res = plus(summand1,summand2) % overloaded + operator
    pZ = polyZonotope(I) % conversion to polyZonotope object
    res = power(base,exponent) % overloaded .^ operator
    res = prod(I,varargin) % overloaded prod-function
    I = project(I,dims) % projection onto subspace
    I = projectHighDim_(I,N,proj) % project to a high-dimensional space
    I = quadMap(varargin) % quadratic map
    r = rad(I) % radius (half of diameter)
    r = radius(I) % radius of enclosing hyperball
    [res,S] = representsa_(I,type,tol,varargin) % comparison to other representations
    res = rdivide(numerator,denominator) % overloaded ./ operator
    I = reshape(I,varargin) % overloaded reshape-function
    res = sin(I) % sine function
    I = sinh(I) % hyperbolic sine function
    varargout = size(I,varargin) % overloaded size-function
    res = split(I,n) % split along one dimension
    I = sqrt(I) % square root
    I = subsasgn(I,S,val) % value assignment
    newObj = subsref(I,S) % read from object
    res = sum(I,varargin) % overloaded sum-function
    res = supremum(I) % read upper limit
    res = tan(I) % tangent function
    I = tanh(I) % hyperbolic tangent function
    res = times(factor1,factor2) % overloaded .* function
    I = transpose(I) % overloaded .' function
    I = uminus(I) % overloaded unary - operator
    I = uplus(I) % overloaded unary + operator
    I = vertcat(varargin) % vertical concantenation
    zB = zonoBundle(I) % conversion to zonoBundle object
    Z = zonotope(I) % conversion to zonotope object
    
    % display functions
    display(I)
end

methods (Static = true)
    I = generateRandom(varargin) % generates random interval
    I = enclosePoints(points) % enclosure of point cloud
    I = empty(n) % instantiates an empty interval
    I = Inf(n) % instantiates a fullspace interval
end

end


% Auxiliary functions -----------------------------------------------------

function [lb,ub] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 2
        throw(CORAerror('CORA:tooManyInputArgs',2));
    end
    
    % assign lower and upper bound
    [lb,ub] = setDefaultValues({[],[]},varargin);

    % set upper bound to value of lower bound if only one value given
    if isnumeric(lb) && ~isempty(lb) && isempty(ub)
        ub = lb;
    end

end

function aux_checkInputArgs(lb,ub,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0
        
        inputArgsCheck({ ...
            {lb, 'att', 'numeric', 'nonnan'}; ...
            {ub, 'att', 'numeric', 'nonnan'}; ...
        })

        if ~isempty(lb) && ~isempty(ub)
            if ~all(size(lb) == size(ub))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Limits are of different dimension.'));
            elseif length(size(lb)) > 2
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Only 1d and 2d intervals are supported.'));
            elseif ~all(all(lb <= ub))
                % check again using tolerance (little bit slower)
                if ~all(all(lb < ub | withinTol(double(lb),double(ub))))
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Lower limit larger than upper limit.'));
                end
            end
        end
        
    end

end

function [lb,ub] = aux_computeProperties(lb,ub)
% if one dimension is [-Inf,-Inf] or [Inf,Inf], the interval is empty
% (cannot be displayed for interval matrices -> throws error)

    % assign correct size if empty interval
    if isempty(lb) && isempty(ub)
        ub = zeros(size(lb));
    end

    if any(any( isinf(lb) & isinf(ub) & (sign(lb) == sign(ub)) ))
        n = size(lb);
        if all(n > 1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Empty interval matrix cannot be instantiated'));
        end
        n(n==1) = 0;
        lb = zeros(n);
        ub = zeros(n);
    end

end

% ------------------------------ END OF CODE ------------------------------
