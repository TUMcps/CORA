classdef matPolytope
% matPolytope class 
%
% Syntax:
%    obj = matPolytope(V)
%
% Inputs:
%    V - cell-array storing the vertices (n x m x N)
%
% Outputs:
%    obj - generated matPolytope object
%
% Example:
%    V(:,:,1) = [1 2; 0 1];
%    V(:,:,2) = [1 3; -1 2];
%    matP = matPolytope(V);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matZonotope

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       21-June-2010
% Last update:   03-April-2023 (MW, remove property dim)
%                02-May-2024 (TL, new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    V % vertices (n x m x N)

    % legacy
    verts
    vertex
end
    
methods
    
    % class constructor
    function obj = matPolytope(V)

        % 1. check input
        if nargin == 0
            V = [];
        elseif nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end

        % 2. check input arguments
        aux_checkInputArgs(V)

        % 3. compute properties
        V = aux_computeProperties(V);

        % 4. assign properties
        obj.V = V;
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
    res = representsa(matP,varargin)

    % legacy ---

    function h = get.verts(matZ)
        CORAwarning("CORA:deprecated","property",'matPolytope.verts','CORA v2024.2.0','Please use matPolytope/numverts() instead.','This change was made to reduce internal maintenance.');
        h = numgens(matZ);
    end

    function matZ = set.verts(matZ,N)
        CORAerror("CORA:noops",'matZ.verts');
    end
    
    function h = get.vertex(matZ)
        CORAwarning("CORA:deprecated","property",'matPolytope.vertex','CORA v2024.2.0','Please use matPolytope.V instead.','This change was made for clarification, as this properties holds all vertices, not just a single vertex.');
        h = numgens(matZ);
    end

    function matZ = set.vertex(matZ,V)
        CORAerror("CORA:noops",'matZ.verts');
    end
    
end

end


% Auxiliary functions -----------------------------------------------------


function aux_checkInputArgs(V)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED

        if iscell(V)
            % legacy
            inputArgsCheck({ ...
                {V, 'att', 'cell'}; ... % cell is legacy
            })

            % input checks are less rigorous here ..

            for i = 1:numel(V)
                % check each generator matrix
                Vi = V{i};
                if ~all(size(V{1},1:2) == size(Vi,1:2))
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        sprintf('Dimension mismatch between given vertices (e.g. first and %i-th vertex).',i)));  
                end
            end

        else
            inputArgsCheck({ ...
                {V, 'att', 'numeric','nonnan'}; ...
            })
        end
        
    end
end

function V = aux_computeProperties(V)
    if iscell(V)
        % legacy, convert to (n x m x N) shape

        % show warning
        CORAwarning('CORA:deprecated','constructor for matPolytope using a','cell of vertices','CORA v2024.2.0','Please use a single numeric matrix with dimensions (n x m x N) instead.','This change was made to improve speed.')

        % correctly initialize empty vertices
        if isempty(V)
             V = zeros(0,0,0);
             return
        end

        % store given vertices
        V_legacy = V;

        % preallocate new generators
        V = zeros([size(V_legacy{1}),numel(V_legacy)]);

        % copy vertices
        for i=1:numel(V_legacy)
            V(:,:,i) = V_legacy{i};
        end
    else
        % correctly initialize empty vertices
        if isempty(V)
             V = zeros([size(V,1:2),1]);
             return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
