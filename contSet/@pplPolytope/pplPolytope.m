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
    C = 1;
    d = 0;
end
    
methods
    %class constructor
    function obj = pplPolytope(C,d)
        obj.dim = length(C(1,:));
        obj.C = C;
        obj.d = d;
    end
         
    %methods in seperate files 
    obj = and(obj,otherObj)
    V = vertices(obj)
    obj = mtimes(matrix,obj)
    obj = plus(summand1,summand2)
    res = isempty(obj)
    bool_val = iscontained(obj,P)
    obj = project(obj,dims)
    Z = parallelotope(varargin)
    [P] = facelift(P,Z)
    
    %display functions
    plot(varargin)
    display(obj)
end
end

%------------- END OF CODE --------------