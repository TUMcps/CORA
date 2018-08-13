classdef conZonotope < zonotope
% conZonotope - object constructor
%
% Syntax:
%       object constructor: obj = conZonotope(cG, A, b)
%
% Inputs:
%    cG - Matrix containing zonotope center and generators cG = [c,G]
%    A - constraint matrix A*ksi = b
%    b - constraint vector A*ksi = b
%
% Outputs:
%    obj - Generated Object
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%    plotZono(cZono);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Dmitry Grebenyuk
% Written:      03-September-2017
%
% Last update:  ---
%               ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    % constrain A*ksi = b; |ksi| <= 1
    % format:       matrix
    A = [];
    
    % format:       column vector
    b = [];
    
    % the value of ksi at vertexes
    % format:       column vector
    ksi = [];
    
    % R = [rho_l, rho_h] (A.3)
    % format:       column vector
    R = [];
    
end
    
methods
    
    % class constructor
    function obj = conZonotope(cG, A, b)
        
        % cleate a zonotope
        obj@zonotope(cG);
        
        % three inputs: a zonotop and constrain
        if nargin == 3
            obj.A = A;
            obj.b = b;
        elseif nargin == 2 || nargin > 3
            error('Wrong syntax. Type "help conZonotope" for more information.')
        end
        
        
    end
    
    
    % methods in seperate files
    [Z] = and(Z1,Z2)
    res = interval(obj);
    res = isempty(obj);
    res = reduce(obj,method,orderG,orderC,varargin)
    res = rescale(obj,varargin)
    res = mptPolytope(obj);
    res = mtimes(factor1,factor2);
    res = intervalMultiplication(obj,I);
    plot( obj, varargin )
    plotFilled( obj, varargin);
    plotZono(obj,varargin);
    res = project(obj,dim);
    display( obj )
    res = vertices(obj);
    res = verticesKsi(obj);
    res = zonotope(obj);
             
end

% prevent unintensional usage of superclass methods
methods(Hidden)
    
    function res = abs(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = appVolume(obj,r)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = box(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = constrSat(obj,C,d)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = containsPoint(obj,p)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = cubicMultiplication(obj,C)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [c_int, deltaVec] = cubicMultiplication_interval(obj,Cint)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = cubicMultiplication_simple(obj,C)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = deleteAligned(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = deleteZeros(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = dim(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = dirPolytopeNew(obj,direction)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = dominantDirections(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zenclose,rotMatrixInv] = encloseMany(Zdummy,Z,direction)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [P] = enclosingPolytope(varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zenclose] = enclosingZonotope(Zfirst,V,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Z = enlarge(Z,factorVec)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [handle] = eventFcn(obj,x,direction)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Z] = exactPlus(Z1,Z2,varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zrem] = filterOut(Zdummy,Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [genLen] = generatorLength(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function G = generators(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Rguard,Rguard_noInt,split,RbeforeInt] = guardProjection(obj,hSpace,A,B,loc,partialIntersectionFlag,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Rguard,tmin,tmax] = guardProjection_09Sep2013(obj,hSpace,A,B,t_hit,tmin,tmax,loc,partialIntersectionFlag,backwards,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Rguard = guardProjection_constFlow(obj,halfspace,A,B,t_hit,tmin,tmax,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [obj] = halfspace(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function result = in(Z1,Z2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [res] = inORH(Z1,Z2,rotMatrixInv)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [res] = inParallelotope(Z1,Z2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Z] = intersection(Z1,Z2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function result = inViaProj(Z1,Z2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = isIntersecting(obj1,obj2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = isIntersectingApprox(obj1,obj2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Z] = minus(minuend,subtrahend)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zquad] = mixedMultiplication(Z1,Z2,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [P] = multiPolytope(firstZ,Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function res = norm(obj, varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Z = or(Z1, Zcell)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [V] = orthVectors(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [R_guard_noInt, R_new] = partialGuardProjection(obj,hSpace,A,B,loc_old,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Rguard = partialGuardProjection_forward(obj,hSpace,A,B,loc_old,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function plotAsText(varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function p = polygon(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [P,comb] = polytope(Z, varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zquad] = quadraticMultiplication(Z,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zquad] = quadraticMultiplication_interval(Z,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zquad] = quadraticMultiplication_parallel(Z,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Zquad] = quadraticMultiplication_zono(Z,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [min,max] = quadraticProgramming_interval(Z,Q)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [qZ] = quadZonotope(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function  r = radius(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [p] = randPoint(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [p] = randPointExtreme(obj)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Z] = rotate(Z,dims,angle)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [Znew] = splitFirstGen(Zdummy,Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Zres = tensorMultiplication(Z,M,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Zres = tensorMultiplication_zono(Z,M,options)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [V] = underapproximate(varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function Z = unify(Z1, Z2)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function [vol] = volume(Z)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function volumeRatio(varargin)
       error('This operation is not implemented for class "conZonotope"!'); 
    end
    
end


end

%------------- END OF CODE -------