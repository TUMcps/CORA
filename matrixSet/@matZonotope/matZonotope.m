classdef matZonotope
% matZonotope class 
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
% See also: intervalhull,  polytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    dim = 1;
    gens = 0;
    center = 0; 
    generator = [];
end
    
methods
    %class constructor
    function obj = matZonotope(input1,input2)
        %one input
        if nargin==0
            matrixCenter = 0;
            matrixGenerator = [];
        elseif nargin==1
            if isa(input1,'zonotope')
                Z=get(input1,'Z');
                %extract center
                c=Z(:,1);
                %extract generator matrix
                G=Z(:,2:end);
                %obtain matrix center
                matrixCenter = vec2mat(c);
                %obtain matrix generators
                for i=1:length(G(1,:))
                    matrixGenerator{i}=vec2mat(G(:,i));
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
        
    %display functions
    plot(varargin)
    display(obj)

end
end

%------------- END OF CODE --------------