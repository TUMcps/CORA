classdef spectraShadow < contSet
% spectraShadow - object constructor for spectraShadow objects
%
% Description:
%    This class represents spectrahedral shadow objects defined as
%    SpS = {Gx+c|A0 + x1*A1 + ... + xm*Am >= 0}
%      (projection representation)
%    or, equivalently,
%    SpS = {y|\exists z, A0 + y1*B1 + ... + ym*Bm + z1*C1 + ... + zl*Cl>=0}
%      (existential sum representation)
%    where the Ai, Bi, Ci are all real, symmetric matrices, and the ">=0"
%    means that the preceding matrix is positive semi-definite.
%
%
% Syntax:
%    SpS = spectraShadow(A)
%    SpS = spectraShadow(A,c)
%    SpS = spectraShadow(A,c,G)
%    SpS = spectraShadow(ESumRep)
%
% Inputs:
%    A - coefficient matrices, i.e., A = [A0 A1 ... Am]
%    c - center vector
%    G - generator matrix
%    ESumRep - Cell array with two elements, containing the existential sum
%       representation of the spectrahedral shadow, i.e.,
%       ESumRep = {[B0 B1 ... Bm] [C1 ... Cl]}
%
% Outputs:
%    obj - generated spectraShadow object
%
% Example: 
%    % 2 dimensional box with radius 3 around point [-1;2]:
%    A0 = eye(3);
%    A1 = [0 1 0;1 0 0;0 0 0];
%    A2 = [0 0 1;0 0 0;1 0 0];
%    SpS = spectraShadow([A0,A1,A2],[-1;2],3*eye(2));
%
% References:
%    [1] T. Netzer. "Spectrahedra and Their Shadows", 2011
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       18-April-2023 
% Last update:   01-August-2023 (AK, restructured arguments)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % spectrahedron coefficient matrices A = [A0,A1,...,Ad]
    % m+1 matrices of dimension k x k => k x k*(m+1)
    A;

    % generator representation of the spectrahedral shadow
    % -> Generator matrix G + center vector c
    % n x m matrix, n x 1 vector
    G;
    c;
end

properties (SetAccess = protected, GetAccess = public)
    % emptiness
    emptySet;

    % full-dimensionality
    fullDim;

    % boundedness
    bounded;
    
    % center of the 'base' spectrahedron
    center;
    
    % Existential Sum representation
    % Alternative representation of a spectrahedron needed for certain
    % operations such as the intersection
    ESumRep;
end

methods

    function obj = spectraShadow(varargin)
        
        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end
        assertNarginConstructor(1:4,nargin);
        
        % 0. init hidden properties
        obj.emptySet = setproperty();
        obj.fullDim = setproperty();
        obj.bounded = setproperty();
        obj.center = setproperty();
        obj.ESumRep = setproperty();
        
        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'spectraShadow')
            SpS = varargin{1};
            obj.A = SpS.A; obj.c = SpS.c; obj.G = SpS.G;
            obj.precedence = SpS.precedence;
            
            obj.ESumRep.val = SpS.ESumRep.val;
            obj.emptySet.val = SpS.emptySet.val;
            obj.fullDim.val = SpS.fullDim.val;
            obj.bounded.val = SpS.bounded.val;
            obj.center.val = SpS.center.val;
            return
        end

        % 2. parse input arguments: varargin -> vars
        [A,c,G,ESumRep] = aux_parseInputArgs(varargin{:});
        % Note that at this stage, A may be a cell containing ESumRep

        % 3. check correctness of input arguments
        aux_checkInputArgs(A,c,G,ESumRep,nargin);
        
        % 4. compute properties and hidden properties
        [A,c,G,ESumRep] = aux_computeProperties(A,c,G,ESumRep,nargin);

        % 5. assign properties
        obj.A = sparse(A);
        obj.c = sparse(c);
        obj.G = sparse(G);
        obj.ESumRep.val = ESumRep;

        % 6. set precedence (fixed)
        obj.precedence = 40;
    end
end

methods (Static = true, Access = public)
    SpS = generateRandom(varargin); % Generate random spectrahedral shadow
    SpS = empty(varargin); % Instantiate empty spectrahedron
    SpS = Inf(varargin); % Instantiate fullspace spectrahedron
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [A,c,G,ESumRep] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables
    
    % Set default values
    [A,c,G] = setDefaultValues({0,[],[]},varargin);
    
    % Identify if initialization is made via just A or the existential sum
    % representation
    if nargin == 1 && iscell(A)
        ESumRep = A;
        A = 0;
    else
        ESumRep = [];
    end

end

function aux_checkInputArgs(A,c,G,ESumRep,n_in)
% check correctness of input arguments

    if CHECKS_ENABLED

            % First, we need to check whether the first argument is a cell; if it
            % is, we just set ESumRep and stop right there. Otherwise, we continue
            % with the 'usual' set up of the spectrahedron
            if n_in == 1 && ~isempty(ESumRep)
                % To construct the corresponding A,c,G, we need to check that
                % ESumRep does indeed consist of matrices

                if CHECKS_ENABLED

                    if (size(ESumRep,1) * size(ESumRep,2)) ~= 2
                        % (That was a little hack to make sure we have a (2,1) or
                        % (1,2) cell array)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The cell array for ESumRep is not a (2,1) or (1,2) array.'));
                    end

                    if ~ismatrix(ESumRep{1}) || ~ismatrix(ESumRep{2})
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'One of the entries of ESumRep is not a matrix.'));
                    end

                    if size(ESumRep{1},1) ~= size(ESumRep{2},1) && ~isempty(ESumRep{1}) && ~isempty(ESumRep{2})
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The dimensions of the two matrices in ESumRep do not match.'));
                    end

                end


                % We can now safely construct A
                A = [ESumRep{1} ESumRep{2}];

                % For c and G, we need to proceed similarly to the 'usual' nargin=1
                % case, the only difference being that we need G = [I 0], instead
                % of G = I

                if CHECKS_ENABLED
                    % First, need to check that A is a matrix
                    inputArgsCheck({ ...
                        {A, 'att', 'numeric', 'nonnan'}; ...
                        })

                end


                % Deduce the dimensions
                k1 = size(ESumRep{1},1);
                m1 = size(ESumRep{1},2) / k1 - 1;
                k2 = size(ESumRep{2},1);
                m2 = size(ESumRep{2},2) / k1;

                if CHECKS_ENABLED

                    % Check that k1 and k2 are equal (or that k2 is zero)
                    if k2~=0 && k1~=k2
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The coefficient matrix ESumRep{1} does not have the same height as ESumRep{2}.'));
                    end

                    % Check that m1 and m2 are integers
                    if (mod(m1,1) ~= 0) && (m1 ~= 0)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The coefficient matrix ESumRep{1} does not have the right dimension (should be k x k*(m+1)).'));
                    end
                    if mod(m2,1) ~= 0 && (m2 ~= 0)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The coefficient matrix ESumRep{2} does not have the right dimension (should be k x k*(m+1)).'));
                    end

                    
                end

                if m2 ~= 0
                    G = [speye(m1) sparse(m1,m2)];
                else
                    G = speye(m1);
                end
                c = zeros([m1 1]);


            elseif n_in == 1 && isempty(ESumRep)

                % If nargin = 1, we need to divide size(A,2) by size(A,1). For
                % that, we need to make sure that this is allowed, i.e., that A is
                % a matrix, and that the dimensions somehow match:

                if CHECKS_ENABLED
                    % First, need to check that A is a matrix
                    inputArgsCheck({ ...
                        {A, 'att', 'numeric', 'nonnan'}; ...
                        })
                end


                % Deduce the dimensions
                k = size(A,1);
                m = size(A,2) / k - 1;

                if CHECKS_ENABLED
                    % Check that m is an integer
                    if (mod(m,1) ~= 0) && (m ~= 0)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The coefficient matrix A does not have the right dimension (should be k x k*(m+1)).'));
                    end
                end

                G = speye(m);
                c = zeros(m,1);
            elseif n_in == 2
                G = speye(size(c, 1));
            end

        if n_in > 0

            % Checking types
            inputArgsCheck({ ...
                {A, 'att', 'numeric', {'nonnan', 'nonempty'}}; ...
                {c, 'att', 'numeric', 'nonnan'}; ...
                {G, 'att', 'numeric', 'nonnan'}; ...
                })
            if ~ismatrix(A)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Coefficient matrix A is not a matrix.'));
            elseif ~ismatrix(G)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Generator matrix G is not a matrix.'));
            elseif ~isvector(c) && ~isempty(c)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Center vector c is not a vector.'));
            end

            % To deal with empty cases, we need a dedicated emptyness-check
            isemptyMatrix = @(M) (size(M,1) == 0) && (size(M,2) == 0);

            % The matrix A is not allowed to be empty, at all
            if isemptyMatrix(A)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The matrix A is empty. If you are trying to create an empty or a fullspace spectrahedral shadow, please consider calling spectraShadow.empty or spectraShadow.Inf instead.'));
            end

            if isemptyMatrix(c) && isemptyMatrix(G) && (size(A,1) == size(A,2))
                % Special case: A is non-empty, but there is not enough
                % data to deduce what dimension the spectrahedron lives in
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Both the generator matrix G and the center vector c ',...
                    'are empty, while the coefficient matrix A is only a ',...
                    'k x k-matrix. It is thus impossible to deduce in ',...
                    'which space the spectrahedron lives in ',...
                    '(i.e., for which n there holds that S is in R^n).']));
            end

            if isemptyMatrix(c) && ~isemptyMatrix(G)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                        'c is empty, but G is not. This is not consistent.'));
            end

            % Deduce the dimensions
            k = size(A,1);
            m = size(A,2) / k - 1;
            n = size(c,1);

            % Check that m is an integer
            if (mod(m,1) ~= 0) && (m ~= 0)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The coefficient matrix A does not have the right dimension (should be k x k*(m+1)).'));
            end

            % Check dimensions
            if size(G,1) ~= n
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Dimension mismatch between the center vector c and the generator matrix G.'));
            elseif size(G,2) ~= m
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Dimension mismatch between the coefficient matrix A and the generator matrix G.'));
            end

            % Check that each submatrix of A is symmetric
            for i = 1:(m+1)
                Ai = A(:,k*(i-1)+1:k*i);
                tol = 1e-6;
                if ~isApproxSymmetric(Ai,tol)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'The coefficient matrix A does not consist of symmetric sub-matrices.'));
                end
            end
        end
    end

end

function [A,c,G,ESumRep] = aux_computeProperties(A,c,G,ESumRep,n_in)
% compute properties
    if n_in == 1 && ~isempty(ESumRep)
        A = [ESumRep{1} ESumRep{2}];
        
        k = size(A,1);
        m1 = size(ESumRep{1},2)/k-1;
        m2 = size(ESumRep{2},2)/k;
        
        G = [speye(m1) sparse(m1,m2)];
        c = sparse(m1,1);
    elseif n_in == 1 && isempty(ESumRep)
        k = size(A,1);
        m = size(A,2)/k-1;
        G = speye(m);
        c = sparse(m, 1);
        ESumRep = {[] []};
    elseif n_in == 2
        G = speye(size(c,1));
        ESumRep = {[] []};
    elseif n_in == 3
        ESumRep = {[] []};
    end
end

% ------------------------------ END OF CODE ------------------------------
