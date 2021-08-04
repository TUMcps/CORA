classdef (InferiorClasses = {?mp}) interval < contSet
% interval - object constructor for real-valued intervals 
%
% Description:
%    This class represents interval objects defined as
%    {x | a_i <= x <= b_i, \forall i = 1,...,n}.
%
% Syntax:
%       obj = interval()
%       obj = interval(I)
%       obj = interval(a)
%       obj = interval(a,b)
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      19-June-2015
% Last update:  18-November-2015
%               26-January-2016
%               15-July-2017 (NK)
%               01-May-2020 (MW, delete redundant if-else)
%               20-March-2021 (MW, errConstructor)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    inf (:,:) double {mustBeNumeric} = [];
    sup (:,:) double {mustBeNumeric} = [];
end

methods
    %class constructor
    function obj = interval(varargin)
        
        % one input argument
        if nargin==1
            
            if isa(varargin{1},'interval')
                % copy constructor
                obj = varargin{1};
            elseif isnumeric(varargin{1})
                obj.inf = varargin{1};
                obj.sup = varargin{1};
            else
                [id,msg] = errConstructor(); error(id,msg);
            end
        
        % two input arguments
        elseif nargin==2
            
            % check sizes
            if any(any(size(varargin{1}) - size(varargin{2})))
                [id,msg] = errConstructor('Limits are of different dimension.');
                error(id,msg);
            elseif all(all(varargin{1} <= varargin{2}))
                obj.inf = varargin{1};
                obj.sup = varargin{2};
            else
                [id,msg] = errConstructor('Lower limit larger than upper limit.');
                error(id,msg);
            end
        
        elseif nargin > 2
            
            % too many input arguments
            [id,msg] = errConstructor('Too many input arguments.'); error(id,msg);
            
        end
        
        % set parent object properties
        obj.dimension = size(obj.inf,1);
    end
    
    function ind = end(obj,k,n)
    % overloads the end operator for referencing elements, e.g. I(end,2),
        ind = size(obj,k);
    end
     
    %methods in seperate files 
    res = plus(summand1,summand2)
    res = minus(minuend,subtrahend)
    res = mtimes(factor1,factor2)
    res = mrdivide(numerator,denominator)
    res = mpower(base,exponent)
    intVal = uplus(intVal)
    intVal = uminus(intVal) 
    res = exp(exponent) %exponential function
    res = abs(value) %absolute value function
    res = sin(intVal) %sine function
    res = cos(intVal) %cosine function
    res = tan(intVal) %tangent function
    newObj = subsref(obj, S) % retrieves values from arrays
    obj = subsasgn(obj, S, value) % assigns values to arrays
    obj = horzcat(varargin) % horizontal concatenation (a = [b,c])
    obj = vertcat(varargin) % vertical concatenation (a = [b;c])
    res = isscalar(intVal) % checks if interval is scalar
    res = supremum(obj) % returns supremum
    res = infimum(obj) % returns infimum
    res = center(obj) % returns center
    res = dim(obj) % return dimension
    res = rad(obj) % returns radius
    obj = and(obj,otherObj) % intersection
    res = length(obj) % returns the length of the array
    res = sqrt(obj) % returns the square root
    varargout = size(obj, varargin) % returns size of object
    res = sinh(intVal) %hyperbolic sine function
    res = cosh(intVal) %hyperbolic cosine function
    res = tanh(intVal) %hyperbolic tangent function
    res = at(intVal, index) % a(i, j)
    res = isIntersecting(obj1,obj2,varargin)
    res = split(input, number) % devides an interval by two
    coordinateMat = gridPoints(obj,segments) % compute uniformly distributed points
    res = in(obj1,obj2) %determines if a set obj2 is entirely in an interval obj1
    res = le(obj1,obj2) %Overloads the <= operator; here: Is one interval equal or the subset of another interval?
    res = lt(obj1,obj2) %Overloads the < operator; here: Is one interval the subset of another interval?
    P = polytope(obj,options) %Converts an interval object to a polytope object
    Z = zonotope(obj) %Converts an interval object into a zonotope object
    r = radius(obj) %Computes radius of enclosing hyperball of an interval 
    V = vertices(obj) %Computes vertices of an interval object
    V = volume(obj) % Computes volume of an interval
    obj = reshape(varargin) %Overloads the opertor 'reshape' for reshaping matrices
    res = sum(intVal) %Overloaded 'sum()' operator for intervals
    h = plot(varargin) %Plots 2-dimensional projection of an interval
    intVal = ctranspose(intVal) %Overloaded ''' operator for single operand
    res = ne( int1, int2 ) % ' ~= ' overloading
    res = prod( obj, dim ) % product of array elements
    res = diag(obj) % create diagonal matrix or get diagonal elements of matrix
    [value,isterminal,direction] = eventFcn(obj,x,direction)
    val = get(obj, propName)
    res = convHull(Int1,Int2) % convex hull of two intervals
    pZ = polyZonotope(obj)
    
    %display functions
    disp(obj)
end

methods (Static = true)
    Int = generateRandom(varargin) % generates random interval
    Int = enclosePoints(points) % enclose point cloud with inteterval
end


end

%------------- END OF CODE -------