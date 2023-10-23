classdef (InferiorClasses = {?mp}) matZonotope
% matZonotope class 
%
% Syntax:
%    obj = matZonotope()
%    obj = matZonotope(C,G)
%
% Inputs:
%    C - center matrix
%    G - cell-array storing the generator matrices
%
% Outputs:
%    obj - generated matZonotope object
%
% Example:
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2];
%    G{2} = [2 0; 1 -1];
%
%    matZ = matZonotope(C,G);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matPolytope

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
%                04-June-2010
%                27-August-2019
%                03-April-2023 (MW, remove property dim)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    gens = 0;
    center = []; 
    generator = [];
end
    
methods
    
    % class constructor
    function obj = matZonotope(input1,input2)

        if nargin == 0
            % empty matrix zonotope
            matrixCenter = [];
            matrixGenerator = [];
            nrGens = 0;

        elseif nargin == 1
            if isa(input1,'zonotope')
                % conversion from zonotope
                % extract center and generator matrix
                c = center(input1);
                G = generators(input1);
                % obtain matrix center
                matrixCenter = c;
                % obtain matrix generators
                if isempty(G)
                    matrixGenerator = [];
                    nrGens = 0;
                elseif isvector(G)
                    nrGens = 1;
                    matrixGenerator = {G};
                else
                    nrGens = size(G,2);
                    matrixGenerator = num2cell(G,[1,nrGens]);
                end
            elseif isa(input1,'matZonotope')
                % copy constructor
                nrGens = input1.gens;
                matrixCenter = input1.center;
                matrixGenerator = input1.generator;
            else
                % only center given...
                matrixCenter = input1;
                matrixGenerator = [];
                nrGens = 0;
            end

        elseif nargin == 2

            % TODO: integrate dimension check...
            matrixCenter = input1;
            matrixGenerator = input2;
            nrGens = length(input2);

        else
            
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end

        % set properties
        obj.gens = nrGens;
        obj.center = matrixCenter;
        obj.generator = matrixGenerator;
    end
         
    %methods in seperate files     
    matZ = plus(summand1,summand2)
    matZ = mtimes(factor1,factor2)
    matZ = mpower(matZ,exponent)
    matZ = powers(matZ,varargin)
    matZ = expmInd(matZ,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmIndMixed(matZ,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,u)
    intMat = intervalMatrix(matZ,varargin)
    matZ = zonotope(matZ)
    dist = expmDist(matZ,intMat,maxOrder)
    matZred = reduce(matZ,option,order,filterLength)
    vol = volume(matZ)
    matZ1 = concatenate(matZ1,matZ2)
    res = norm(matZ,varargin)
    newObj = subsref(matZ,S)
        
    %display functions
    plot(matZ,varargin)
    display(matZ)

end
end

% ------------------------------ END OF CODE ------------------------------
