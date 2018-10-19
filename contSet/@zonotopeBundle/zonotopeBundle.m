classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotopeBundle
%classdef(InferiorClasses = {?intervalMatrix, ?matZonotope, ?matPolytope, ?intval}) zonotopeBundle
% zonotopeBundle class 
%
% Syntax:  
%    object constructor: Obj = zonotopeBundle(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - cell array of zonotopes
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
% Written:      09-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    Z = [];
    parallelSets = 0;
end
    
methods
    %class constructor
    function obj = zonotopeBundle(input)
        
        %one input
        if nargin==1
            %set zonotope cell array
            obj.Z=input;
            %get number of parallel sets
            obj.parallelSets = length(input);
        end
    end
         
    %methods in seperate files     
    Zbundle = plus(summand1,summand2)
    Zbundle = mtimes(factor1,factor2)
    IH = interval(Zbundle)
    Zbundle = reduce(Zbundle,varargin)
    Zred = reduceCombined(Zbundle,option,varargin)
    Zbundle1=enclose(Zbundle1,Zbundle2)
    Zbundle=and(Zbundle1,Zbundle2)
    [P] = enclosingPolytope(varargin)
    [P] = polytope(varargin)
    [P] = parallelotope(varargin)
    res = conZonotope(obj)
    Z1 = cartesianProduct(Z1,Z2)
    Zbundle = enlarge(Zbundle,factorVec)
    Zbundle = shrink(Zbundle,filterLength)
    Zbundle = shrink2(Zbundle,W)
    Zbundle = shrinkIH(Zbundle)
    Zbundle = replace(Zbundle,index,Z)
    Zbundle = shrink3(Zbundle,W,filterLength)
    Zbundle = encloseTight(Zbundle1,Zbundle2,W)
    [vol] = volume(Zbundle)
    [c] = center(Z)
    [Zsplit] = split(Zbundle,options,varargin)
    [Zbundle] = project(Zbundle,dim)
    [Rguard,Rguard_noInt,split,RbeforeInt] = guardProjection(obj,hSpace,A,B,loc,partialIntersectionFlag,options)
    [Rguard] = partialGuardProjection(obj,halfspace,A,B,loc_old,options)
    [Z] = quadraticMultiplication(Z,Q)
    [Zbundle] = quadraticMultiplication_zono(Zbundle,Q)
        
    %display functions
    handle = plot(varargin)
    display(obj)

end
end

%------------- END OF CODE --------------