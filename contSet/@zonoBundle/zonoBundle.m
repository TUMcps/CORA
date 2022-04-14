classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonoBundle  < contSet
% zonoBundle class 
%
% Syntax:  
%    obj = zonoBundle(list)
%
% Inputs:
%    list - cell-array list = {Z1,...,Z2} storing the zonotopes that
%           define the zonotope bundle
%
% Outputs:
%    obj - generated object
%
% Example:
%    zono1 = zonotope([1 3 0; 1 0 2]);
%    zono2 = zonotope([0 2 2; 0 2 -2]);
%
%    zB = zonoBundle({zono1,zono2});
%
%    figure
%    hold on
%    plot(zB,[1,2],'r','Filled',true);
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'g');
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    Z (:,1) cell;
    parallelSets {mustBeNonnegative} = 0;
end
    
methods
    %class constructor
    function obj = zonoBundle(varargin)
        
        %one input
        if nargin==1
            if isa(varargin{1},'zonoBundle')
                % copy constructor
                obj = varargin{1};
            else
                %set zonotope cell array
                obj.Z=varargin{1};
                %get number of parallel sets
                obj.parallelSets = length(varargin{1});
            end
        end
        
        % set parent object properties
        if isempty(obj.Z)
            obj.dimension = 0;
        else
            obj.dimension = size(obj.Z{1}.Z,1);
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
    Z1 = cartProd(Z1,Z2)
    Zbundle = enlarge(Zbundle,factorVec)
    Zbundle = shrink(Zbundle,filterLength)
    Zbundle = replace(Zbundle,index,Z)
    Zbundle = encloseTight(Zbundle1,Zbundle2,W)
    [vol] = volume(Zbundle)
    [c] = center(Z)
    [Zsplit] = split(Zbundle,options,varargin)
    [Zbundle] = project(Zbundle,dim)
    [Z] = quadMap(Z,Q)
    [Zbundle] = quadMap_zono(Zbundle,Q)
        
    %display functions
    handle = plot(varargin)
    display(obj)

end

methods (Static = true)
    Zbundle = generateRandom(varargin) % generate random zonotope bundle
end


end

%------------- END OF CODE --------------