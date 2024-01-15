classdef (InferiorClasses = {?mp}) intervalMatrix 
% intervalMatrix class 
%
% Syntax:
%    obj = intervalMatrix()
%    obj = intervalMatrix(C,D)
%
% Inputs:
%    C - center matrix
%    D - width matrix
%
% Outputs:
%    obj - generated object
%
% Example:
%   C = [0 2; 3 1];
%   D = [1 2; 1 1];
%   intMat = intervalMatrix(C,D);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       18-June-2010
% Last update:   26-August-2011
%                15-June-2016
%                06-May-2021
%                03-April-2023 (MW, remove properties dim and setting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % interval
    int = [];
end
    
methods 
    %class constructor
    function obj = intervalMatrix(matrixCenter,matrixDelta)

        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));

        elseif nargin == 1
            
            if isa(matrixCenter,'intervalMatrix')
                % copy constructor
                obj = matrixCenter;
            else
                % only center given, radius = 0
                obj.int = interval(matrixCenter,matrixCenter);
            end

        elseif nargin==2
            % ensure positive matrix deltas
            matrixDelta=abs(matrixDelta);
            obj.int = interval(matrixCenter-matrixDelta,matrixCenter+matrixDelta);
        end
    end
         
    %methods in seperate files 
    intMat = plus(summand1,summand2)
    intMat = mtimes(factor1,factor2)
    intMat = mpower(intMat,exponent)
    intMat = powers(intMat,varargin)
    [eI, iPow, E] = expm(intMat,varargin)
    [eI, iPow, E] = expmInd(intMat,varargin)
    [eI,eI2,iPow,iPow2,E] = expmMixed(intMat,r,intermediateOrder,maxOrder)
    [eI,eI2,iPow,iPow2,E] = expmIndMixed(intMat,intermediateOrder,maxOrder)
    M = abs(intMat)
    n = infNorm(intMat)
    E = exponentialRemainder(intMat,maxOrder)
    I = interval(intMat)
    V = vertices(intMat)
    matP = matPolytope(intMat)
    matZ = matZonotope(intMat)
    dist = expmDist(intMat,exactMat,maxOrder)
    vol = volume(intMat)
    val = expmNorm(intMat,t)
    val = expmNormInf(intMat,t)
    absBound = expmAbsoluteBound(intMat,t)
    normBoundErr = expmNormErr(intMat,r)  
    normBoundErr = expmNormErrInf(intMat,r)
    element = subsref(intMat, S)
    sq = exactSquare(intMat)
    res = norm(intMat,varargin)
    res = isempty(intMat)
    M = randPoint(intMat,varargin)
    res = contains(intMat,M,type,tol)
    
    %display functions
    plot(varargin)
    display(intMat)
end

methods (Static = true)
    intMat = generateRandom(varargin) % generates random interval matrix
end

end

% ------------------------------ END OF CODE ------------------------------
