classdef (InferiorClasses = {?mp}) matZonotope
% matZonotope class 
%
% Syntax:  
%   obj = matZonotope(C,G)
%
% Inputs:
%    C - center matrix
%    G - cell-array storing the generator matrices
%
% Outputs:
%    obj - generated object
%
% Example:
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2];
%    G{2} = [2 0; 1 -1];
%
%    mz = matZonotope(C,G);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matPolytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
%               27-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    dim = 1;
    gens = 0;
    center = 0; 
    generator = [];
end
    
methods
    
    % class constructor
    function obj = matZonotope(input1,input2)
        %one input
        if nargin==0
            matrixCenter = 0;
            matrixGenerator = [];
        elseif nargin==1
            if isa(input1,'zonotope')
                %extract center
                c=center(input1);
                %extract generator matrix
                G=generators(input1);
                %obtain matrix center
                matrixCenter = vec2mat(c);
                %obtain matrix generators
                if ~isempty(G)
                    for i=1:length(G(1,:))
                        matrixGenerator{i}=vec2mat(G(:,i));
                    end
                else
                    matrixGenerator{1} = zeros(size(matrixCenter));
                end
            else
                matrixCenter=input1;
                matrixGenerator = [];
            end
        elseif nargin==2
            matrixCenter = input1;
            matrixGenerator = input2;
        end
        %set parameters
        obj.dim = length(matrixCenter);
        obj.gens = length(matrixGenerator);
        obj.center = matrixCenter;
        obj.generator = matrixGenerator;
    end
         
    %methods in seperate files     
    matZ = plus(summand1,summand2)
    matZ = mtimes(factor1,factor2)
    matZ = mpower(matZ,exponent)
    matZ = powers(varargin)
    matZ = expmInd(matZ,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmIndMixed(matZ,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,u)
    intMat = intervalMatrix(varargin)
    matZ = zonotope(matZ)
    dist = expmDist(matZ,intMat,maxOrder)
    matZred = reduce(matZ,option,order,filterLength)
    vol = volume(matI)
    matZ1 = concatenate(matZ1,matZ2)
    res = norm(obj, varargin)
    newObj = subsref(obj, S)
        
    %display functions
    plot(varargin)
    display(obj)

end
end

%------------- END OF CODE --------------