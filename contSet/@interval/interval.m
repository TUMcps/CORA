classdef interval
% interval class 
%
% Syntax:  
%    object constructor: Obj = interval(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - left limit
%    input2 - right limit
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      19-June-2015
% Last update:  18-November-2015
%               26-January-2016
%               15-July-2017 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    inf = [];
    sup = [];
    setting = 'sharpivmult';
end
    
methods 
    %class constructor
    function obj = interval(leftLimit,rightLimit)
        
        %no input
        if nargin==0
            obj.inf = [];
            obj.sup = [];
        
        %one input
        elseif nargin==1
            if isa(leftLimit,'interval')
                obj = leftLimit;
            elseif isnumeric(leftLimit)
                obj.inf = leftLimit;
                obj.sup = leftLimit;
            end
        %two inputs
        elseif nargin==2
            if isscalar(leftLimit)
                if leftLimit <= rightLimit
                    obj.inf = leftLimit;
                    obj.sup = rightLimit;
                else
                    disp('Left limit larger than right limit');
                end
            else
                if all(all(leftLimit <= rightLimit))
                    obj.inf = leftLimit;
                    obj.sup = rightLimit;
                else
                    disp('Left limit larger than right limit');
                end
            end
        end
%         % Register the variable as an object
%         obj = class(obj, 'interval');
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
    res = mid(obj) % returns center
    res = rad(obj) % returns radius
    obj = and(obj,otherObj) % intersection
    res = length(obj) % returns the length of the array
    res = sqrt(obj) % returns the square root
    varargout = size(obj, varargin) % returns size of object
    res = sinh(intVal) %hyperbolic sine function
    res = cosh(intVal) %hyperbolic cosine function
    res = tanh(intVal) %hyperbolic tangent function
    res = at(intVal, index) % a(i, j)
    obj = hull(intVal1, intVal2) %the union of two intervals
    res = intersect(int1,int2) %the intersection of two intervals
    res = isIntersectingApprox(obj1,obj2)
    res = isIntersecting(obj1,obj2)
    res = split(input, number) % devides an interval by two
    coordinateMat = gridPoints(obj,segments) % compute uniformly distributed points
    res = in(obj1,obj2) %determines if a set obj2 is entirely in an interval obj1
    res = le(obj1,obj2) %Overloads the <= operator; here: Is one interval equal or the subset of another interval?
    res = lt(obj1,obj2) %Overloads the < operator; here: Is one interval the subset of another interval?
    P = polytope(obj,options) %Converts an interval object to a polytope object
    Z = zonotope(obj) %Converts an interval object into a zonotope object
    r = enclosingRadius(obj) %Computes radius of enclosing hyperball of an interval 
    V = vertices(obj) %Computes vertices of an interval object
    V = volume(obj) % Computes volume of an interval
    obj = reshape(varargin) %Overloads the opertor 'reshape' for reshaping matrices
    res = sum(intVal) %Overloaded 'sum()' operator for intervals
    plotFilled(varargin) %Plots 2-dimensional projection of an interval which is colored inside
    intVal = ctranspose(intVal) %Overloaded ''' operator for single operand
    res = ne( int1, int2 ) % ' ~= ' overloading
    res = prod( obj, dim ) % product of array elements
    res = diag(obj) % create diagonal matrix or get diagonal elements of matrix
    [value,isterminal,direction] = eventFcn(obj,x,direction)
    val = get(obj, propName)
    
    %display functions
    disp(obj)
end
end

%------------- END OF CODE -------