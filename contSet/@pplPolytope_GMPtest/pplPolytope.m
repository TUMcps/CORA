classdef pplPolytope
% pplPolytope class 
%
% Syntax:  
%    object constructor: obj = pplPolytope(varargin)
%    copy constructor: obj = otherObj
%
% Inputs:
%    input1 - C matrix
%    input2 - d vector
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    dim = 1;
    precision = 30;
    C = mp(1,precision);
    d = mp(0,precision);
end
    
methods
    %class constructor
    function obj = pplPolytope(C,d,precision)
        %one input
        if nargin == 2
            obj.dim = length(C(1,:));
            obj.C =  mp(C,precision);
            obj.d = mp(d,precision);
        elseif nargin == 3
            obj.dim = length(C(1,:));
            obj.precision = precision;
            obj.C =  mp(C,precision);
            obj.d = mp(d,precision);
        end
    end
         
    %methods in seperate files 
    obj = and(obj,otherObj)
    V = vertices(obj)
    obj = mtimes(matrix,obj)
    obj = plus(summand1,summand2)
    res = isempty(obj)
    obj = project(obj,dims)
    
    %display functions
    plot(varargin)
    display(obj)
end
end

%------------- END OF CODE --------------