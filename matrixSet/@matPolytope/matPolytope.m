classdef matPolytope
% matPolytope class 
%
% Syntax:
%    obj = matPolytope(V)
%
% Inputs:
%    V - cell-array storing the vertices
%
% Outputs:
%    obj - generated matPolytope object
%
% Example:
%    V{1} = [1 2; 0 1];
%    V{2} = [1 3; -1 2];
%    matP = matPolytope(V);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matZonotope

% Authors:       Matthias Althoff
% Written:       21-June-2010
% Last update:   03-April-2023 (MW, remove property dim)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    verts = 0; % number of vertices
    vertex = [];
end
    
methods
    
    % class constructor
    function obj = matPolytope(input)
        if nargin == 0
            % take default values

        elseif nargin == 1
            if isa(input,'polytope')
                %get vertices from polytope class
                V=vertices(input);
                %rewrite vertices as matrices
                for i=1:length(V(:,1))
                    matrixVertex{i} = vec2mat(V(i,:));
                end
            else
                matrixVertex = input;
            end
            obj.verts = length(matrixVertex);
            obj.vertex = matrixVertex;

        else
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end
    end
         
    % methods in seperate files
    display(matP) % display on command window
    matP = expmInd(matP,maxOrder)
    [eP,eI] = expmIndMixed(matP,intermediateOrder,maxOrder)
    val = expmDist(matP,exactMat,maxOrder) % deprecated?
    intMat = intervalMatrix(matP) % conversion to interval matrix
    matZ = matZonotope(matP) % conversion to matrix polytope
    matP = mpower(matP,exponent) % exponentiation
    matP = polytope(matP) % conversion to polytope
    matP = mtimes(factor1,factor2) % linear map
    plot(matP,varargin) % plot
    matP = plus(summand1,summand2) % Minkowski addition
    matP = simplePlus(summand1,summand2) % ?
    [r,c] = size(matP,rc) % read out dimension of vertices
    matV = vertices(matP) % read vertices of matrix polytope
    res = isempty(matP) % emptiness check
    [r,c] = dim(matP,rc) % dimension of matrix polytope
    
end

end

% ------------------------------ END OF CODE ------------------------------
