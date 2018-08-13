classdef constrainedHyperplane
% constrainedHyperplane class 
%
% Syntax:  
%    object constructor: Obj = constrainedHyperplane(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    h - cell array of halfspaces
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    h = [];
    C = [];
    d = [];
end
    
methods
    %class constructor
    function obj = constrainedHyperplane(input1, input2, input3)
        %input is a hyperplane and a constraint system Cx<=d
        obj.h = input1;
        obj.C = input2;
        obj.d = input3;
    end
         
    %methods in seperate files     
    [value,isterminal,direction] = eventFcn(obj,x,direction)
    val = get(obj, propName)
    res = zonoIntersect(obj, Z)
    res = zonoIn(obj, Z)
    Rguard = guardProjection(obj,A,B,t_hit,tmin,tmax,R0,options)
    res = isIntersecting(obj1,obj2)
    [result] = in(obj1,obj2)
    [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
        
    %display functions
    display(obj)

end
end

%------------- END OF CODE --------------
