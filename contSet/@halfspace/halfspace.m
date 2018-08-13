classdef halfspace
% halfspace class 
%
% Syntax:  
%    object constructor: Obj = halfspace(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    c - normal vector of the halfspace
%    d - distance to the origin
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
% Written:      06-June-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    c = 1;
    d = 0;
end
    
methods
    %class constructor
    function obj = halfspace(input1,input2)
        if nargin==1
            if isa(input1,'halfspace') %dummy constructor when halfspace function is called
                obj = input1;     
            end
        %one input
        elseif nargin==2
            %set zonotope cell array
            obj.c = input1;
            %get number of parallel sets
            obj.d = input2;
            
            % convert column vector to row vector
            if size(obj.c,1) == 1 
               obj.c = obj.c'; 
            end
        end
    end
         
    %methods in seperate files     
    [value,isterminal,direction] = eventFcn(obj,x,direction)
    val = get(obj, propName)
    [k,L,Q,C] = maps(obj, A, f, x0, R0)
    [k,L,Q,C_c,C_g] = mapsCubic(obj, A, f, x0, R0)
    [k,L,Q,lambdaQuad_zono,Lambda_int] = mapsQuad(obj, A, f, x0, R0, c_tilde, d_tilde)
    [k,L,Q] = mapsQuad_old(obj, A, f, x0, R0)
    res = zonoIntersect(obj, Z)
    [h] = mtimes(factor,obj)
    res = commonPoint(obj, h_other)
    [h] = plus(summand1,summand2)
    res = zonoIn(obj, Z)
    res = zonoPartiallyIn(obj, Z)
    S = project(h, S)
    [tmin,tmax,t_total,Rlast] = intersectionTime_oneDirection(obj,R0,contDynamics,options)
    res = isIntersecting(obj1,obj2)
        
    %display functions
    display(obj)
    plot(obj, dimensions, type)
end
end

%------------- END OF CODE --------------