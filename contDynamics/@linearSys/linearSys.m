classdef linearSys < contDynamics
% linearSys - object constructor for linear time-invariant systems
%
% Description:
%    Generates a linear system object according to the following
%    first-order differential equations:
%       x'(t) = A x(t) + B u(t) + c + E w(t)
%       y(t)  = C x(t) + D u(t) + k + F v(t)
%
% Syntax:
%    obj = linearSys()
%    obj = linearSys(A)
%    obj = linearSys(A,B)
%    obj = linearSys(A,B,c)
%    obj = linearSys(A,B,c,C)
%    obj = linearSys(A,B,c,C,D)
%    obj = linearSys(A,B,c,C,D,k)
%    obj = linearSys(A,B,c,C,D,k,E)
%    obj = linearSys(A,B,c,C,D,k,E,F)
%    obj = linearSys(name,A,B)
%    obj = linearSys(name,A,B,c)
%    obj = linearSys(name,A,B,c,C)
%    obj = linearSys(name,A,B,c,C,D)
%    obj = linearSys(name,A,B,c,C,D,k)
%    obj = linearSys(name,A,B,c,C,D,k,E)
%    obj = linearSys(name,A,B,c,C,D,k,E,F)
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    c - constant input
%    C - output matrix
%    D - feedthrough matrix
%    k - output offset
%    E - disturbance matrix
%    F - output disturbance matrix
%
% Outputs:
%    obj - generated linearSys object
%
% Example:
%    A = [-2 0; 1 -3];
%    B = [1; 1];
%    C = [1 0];
%
%    sys = linearSys(A,B,[],C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       23-January-2007 
% Last update:   30-April-2007
%                04-August-2016 (changed to new OO format)
%                01-November-2017 (constant input added)
%                20-March-2018 (NK, output equation parameter added)
%                07-November-2018 (MA, default values for B and C changed)
%                04-March-2019 (MA, default IDs for inputs and outputs)
%                20-May-2020 (NK, name now optional)
%                19-November-2021 (MW, default values for c, D, k)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                15-January-2023 (MW, allow 0 and 1 input arguments)
%                30-August-2024 (MW, integrate E and F matrices)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A     % system matrix: n x n
    B     % input matrix: n x m
    c     % constant input: n x 1
    C     % output matrix: q x n
    D     % feedthrough matrix: q x m
    k     % output offset: q x 1
    E     % disturbance matrix: n x r
    F     % noise matrix: q x s

    % internally set properties
    taylor = [];    % struct storing values from Taylor expansion
    krylov = [];    % struct storing values for Krylov subspace
end

methods
    
    % class constructor
    function linsys = linearSys(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:9,nargin);
        
        % 1. copy constructor: not allowed due to obj@contDynamics below
        % if nargin == 1 && isa(varargin{1},'linearSys')
        %     obj = varargin{1}; return
        % end

        % 2. parse input arguments: varargin -> vars
        [name,A,B,c,C,D,k,E,F] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,A,B,c,C,D,k,E,F,nargin);

        % 4. compute number of states, inputs, and outputs
        [name,A,B,c,C,D,k,E,F,states,inputs,outputs,dists,noises] = ...
            aux_computeProperties(name,A,B,c,C,D,k,E,F);
        
        % 5. instantiate parent class, assign properties
        linsys@contDynamics(name,states,inputs,outputs,dists,noises);
        linsys.A = A; linsys.B = B; linsys.c = c;
        linsys.C = C; linsys.D = D; linsys.k = k;
        linsys.E = E; linsys.F = F;
%         obj.taylor = taylorLinSys(A);
    end

    % wrapper function to read out auxiliary values stored in sys.taylor,
    % if the requested value is not there, we compute it
    function val = getTaylor(linsys,name,varargin)
        val = computeField(linsys.taylor,name,varargin{:});
    end
end

methods (Static = true)
    linsys = generateRandom(varargin) % generates random linear system
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,c,C,D,k,E,F] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % default name
    def_name = 'linearSys';
    A = []; B = []; c = []; C = []; D = []; k = []; E = []; F = [];

    % parse depending on whether first input argument is the name
    if nargin > 0 && ischar(varargin{1})
        % first input argument: name
        [name,A,B,c,C,D,k,E,F] = setDefaultValues({def_name,A,B,c,C,D,k,E,F},varargin);

    else
        % set default name
        name = def_name;
        [A,B,c,C,D,k,E,F] = setDefaultValues({A,B,c,C,D,k,E,F},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,c,C,D,k,E,F,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % ensure that values have correct data type
        nameCheck = {};
        if ~strcmp(name,'linearSys')
            nameCheck = {{name, 'att', {'char','string'}}};
        end

        inputArgsCheck({ ...
            nameCheck{:}
            {A, 'att', 'numeric', 'square'}
            {B, 'att', 'numeric', 'matrix'}
            {c, 'att', 'numeric'}
            {C, 'att', 'numeric', 'matrix'}
            {D, 'att', 'numeric', 'matrix'}
            {k, 'att', 'numeric'}
            {E, 'att', 'numeric', 'matrix'}
            {F, 'att', 'numeric', 'matrix'}
        });

        % offsets must be vectors
        if ~isempty(c) && ~isvector(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Offset c must be a vector.'));
        elseif ~isempty(k) && ~isvector(k)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Offset k must be a vector'));
        end
        % check dimension of offset for state
        if ~isempty(c) && size(A,1) ~= length(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Length of offset c must match row/column dimension of state matrix A.'));
        end

        % check if input matrix/feedthrough have the correct size
        if ~isempty(B)
            if ~isscalar(B)
                if size(A,1) ~= size(B,1)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Column dimension of input matrix B must match '...
                    'row/column dimension of state matrix A.']));
                end
                inputs = size(B,2);
            else % isscalar(B)
                inputs = size(A,1);
            end
            if ~isempty(D) && ~isscalar(D) && inputs ~= size(D,2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Column dimension of input matrix B must match '...
                    'column dimension of feedthrough matrix D.']));
            end
        end

        % check if disturbance matrix has the correct size
        if ~isempty(E)
            if ~isscalar(E) && size(A,1) ~= size(E,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Column dimension of disturbance matrix E must match '...
                    'row/column dimension of state matrix A.']));
            end
        end

        % check if output matrix/feedthrough matrix have the correct size
        outputs = size(A,1); % assume outputs = states
        if ~isempty(C)
            if ~isscalar(C)
                if size(A,1) ~= size(C,2)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Column dimension of output matrix C must match '...
                        'row/column dimension of state matrix A.']));
                end
                outputs = size(C,1);
            end
            if ~isempty(D) && ~isscalar(D) && outputs ~= size(D,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Row dimension of feedthrough matrix D must match '...
                    'row dimension of output matrix C.']));
            end
            if ~isempty(k) && outputs ~= length(k)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Length of offset k must match row dimension of output matrix C.'));
            end
        end

        % check if noise matrix has the correct size
        if ~isempty(F)
            if ~isscalar(F) && outputs ~= size(F,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Column dimension of noise matrix F must match '...
                    'row/column dimension of state matrix A.']));
            end
        end
        
    end

end

function [name,A,B,c,C,D,k,E,F,states,inputs,outputs,dists,noises] = ...
    aux_computeProperties(name,A,B,c,C,D,k,E,F)
% assign zero vectors/matrices for [] values (from user or by default)
% compute number of states, inputs, outputs, dists, and noises
% difficulty: MATLAB can use 1 instead of eye(n) for multiplication

    % number of states
    states = size(A,1);

    % input matrix and number of inputs
    if ~isempty(A) && isempty(B)
        B = zeros(states,1);
    end
    inputs = states;
    if ~isscalar(B)
        inputs = size(B,2);
    end

    % constant offset
    if isempty(c)
        c = zeros(states,1);
    end

    % output matrix and number of outputs
    if isempty(C)
        C = 1;
    end
    outputs = states;
    if ~isscalar(C)
        outputs = size(C,1);
    end

    % feedthrough matrix
    if isempty(D)
        D = zeros(outputs,inputs);
    end

    % output offset
    if isempty(k)
        k = zeros(outputs,1);
    end

    % disturbance matrix and number of disturbances
    if isempty(E)
        E = zeros(states,1);
    end
    if ~isscalar(E)
        dists = size(E,2);
    else
        dists = states;
    end

    % noise matrix and number of noises
    if isempty(F)
        F = zeros(outputs,1);
    end
    if ~isscalar(F)
        noises = size(F,2);
    else
        noises = outputs;
    end

end

% ------------------------------ END OF CODE ------------------------------
