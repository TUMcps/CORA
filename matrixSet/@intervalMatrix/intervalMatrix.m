classdef intervalMatrix 
% intervalMatrix class 
%
% Syntax:  
%    object constructor: Obj = zonotope(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - zonotope matrix
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      18-June-2010
% Last update:  26-August-2011
%               15-June-2016
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    dim = 1;
    Inf = [];
    Sup = [];
    int = [];
    setting = 'sharpivmult';
end
    
methods 
    %class constructor
    function obj = intervalMatrix(matrixCenter,matrixDelta,setting)

        %one input
        if nargin==1
            if isa(matrixCenter,'intervalMatrix')
                obj = matrixCenter;
            else
                obj.dim = length(matrixCenter);
                obj.Inf = matrixCenter;
                obj.Sup = matrixCenter;
                obj.int = interval(matrixCenter,matrixCenter);
                obj.setting = [];
            end
        elseif nargin==2
            obj.dim = length(matrixCenter);
            %ensure positive matrix deltas
            matrixDelta=abs(matrixDelta);
            obj.Inf = matrixCenter-matrixDelta;
            obj.Sup = matrixCenter+matrixDelta;
            obj.int=interval(obj.Inf,obj.Sup);
            obj.setting = [];
        elseif nargin==3
            obj.dim = length(matrixCenter);
            %ensure positive matrix deltas
            matrixDelta=abs(matrixDelta);
            obj.Inf = matrixCenter-matrixDelta;
            obj.Sup = matrixCenter+matrixDelta;
            obj.int=interval(obj.Inf,obj.Sup);
            obj.setting = setting;  %settings: 'sharpivmult','fastivmult'
        end
    end
         
    %methods in seperate files 
    intMat = plus(summand1,summand2)
    intMat = mtimes(factor1,factor2)
    intMat = mpower(intMat,exponent)
    intMat = powers(varargin)
    [eI, iPow, E] = expm(varargin)
    [eI, iPow, E] = expmInd(varargin)
    [eI,eI2,iPow,iPow2,E] = expmMixed(matI,r,intermediateOrder,maxOrder)
    [eI,eI2,iPow,iPow2,E] = expmIndMixed(matI,intermediateOrder,maxOrder)
    M = abs(intMat)
    n = infNorm(intMat)
    E = exponentialRemainder(intMat,maxOrder)
    IH = interval(intMat)
    V = vertices(intMat)
    matV = dominantVertices(matI, maxNumber)
    matP = matPolytope(intMat)
    matZ = matZonotope(intMat)
    A = randomSampling(intMat,varargin)
    dist = expmDist(intMat,exactMat,maxOrder)
    vol = volume(matI)
    val = expmNorm(intMat,t)
    val = expmNormInf(intMat,t)
    absBound = expmAbsoluteBound(intMat,t)
    normBoundErr = expmNormErr(intMat,r)  
    normBoundErr = expmNormErrInf(intMat,r)
    eI = expmVertex(intMat)
    element = subsref(intMat, S)
    sq = exactSquare(A)
    res = norm(obj, varargin)
    
    %display functions
    plot(varargin)
    display(obj)
end
end

%------------- END OF CODE -------