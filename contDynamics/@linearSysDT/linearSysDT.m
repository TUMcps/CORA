classdef linearSysDT < contDynamics
% linearSysDT - object constructor for linear discrete-time systems
%    
% Description:
%    Generates a discrete-time linear system object according to the 
%    following first-order difference equations:
%       x(k+1) = A x(k) + B u(k) + c + w(k)
%       y(k)   = C x(k) + D u(k) + k + v(k)
%
% Syntax:
%    obj = linearSysDT(A,B,dt)
%    obj = linearSysDT(A,B,c,dt)
%    obj = linearSysDT(A,B,c,C,dt)
%    obj = linearSysDT(A,B,c,C,D,dt)
%    obj = linearSysDT(A,B,c,C,D,k,dt)
%    obj = linearSysDT(name,A,B,dt)
%    obj = linearSysDT(name,A,B,c,dt)
%    obj = linearSysDT(name,A,B,c,C,dt)
%    obj = linearSysDT(name,A,B,c,C,D,dt)
%    obj = linearSysDT(name,A,B,c,C,D,k,dt)
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    c - constant input
%    C - output matrix
%    D - feedthrough matrix
%    k - output offset
%    dt - sampling time
%
% Outputs:
%    obj - generated linearSysDT object
%
% Example:
%    A = [-0.4 0.6; 0.6 -0.4];
%    B = [0; 1];
%    c = [0;0]
%    C = [1 0];
%    dt = 0.4;
%    sys = linearSysDT(A,B,c,C,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       20-March-2020 
% Last update:   14-June-2021 (MA, invoke observe from superclass)
%                19-November-2021 (MW, default values for C, D, k)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A   % system matrix: n x n
    B   % input matrix: n x m
    c   % constant input: n x 1
    C   % output matrix: q x n
    D   % feedthrough matrix: q x m
    k   % output offset: q x 1
    dt  % sampling time
end

methods
    
    % class constructor
    function obj = linearSysDT(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'linearSysDT')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,A,B,c,C,D,k,dt] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,A,B,c,C,D,k,dt,nargin);

        % 4. compute number of states, inputs, and outputs
        [name,A,B,c,C,D,k,dt,states,inputs,outputs] = ...
            aux_computeProperties(name,A,B,c,C,D,k,dt);
        
        % 5. instantiate parent class, assign properties
        obj@contDynamics(name,states,inputs,outputs); 
        obj.A = A; obj.B = B; obj.c = c;
        obj.C = C; obj.D = D; obj.k = k;
        obj.dt = dt;
    end
    
    % invoke function observe so that the superclass can access private
    % functions
    function [R, tcomp] = observe(obj,params,options)
        [R, tcomp] = observe@contDynamics(obj,params,options);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,c,C,D,k,dt] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 8
        throw(CORAerror('CORA:tooManyInputArgs',8));
    end

    % default name
    def_name = 'linearSysDT';
    % init properties
    A = []; B = []; c = []; C = []; D = []; k = []; dt = 0;

    % no input arguments
    if nargin == 0
        name = def_name;
        return
    end

    % last input argument is always sampling time
    dt = varargin{end};
    varargin = varargin(1:end-1);

    % parse depending on whether first input argument is the name
    if ischar(varargin{1})
        % first input argument: name
        [name,A,B,c,C,D,k] = setDefaultValues({def_name,A,B,c,C,D,k},varargin);

    else
        % set default name
        name = def_name;
        [A,B,c,C,D,k] = setDefaultValues({A,B,c,C,D,k},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,c,C,D,k,dt,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % ensure that values have correct data type
        if strcmp(name,'linearSysDT')
            % default name (unless explicitly chosen by user, we have A as
            % first input argument)

            inputArgsCheck({ ...
                {A, 'att', 'numeric', 'square'}
                {B, 'att', 'numeric', 'matrix'}
                {c, 'att', 'numeric'}
                {C, 'att', 'numeric', 'matrix'}
                {D, 'att', 'numeric', 'matrix'}
                {k, 'att', 'numeric'}
                {dt, 'att', 'numeric', {'scalar','positive'}}
            });
            
        else

            inputArgsCheck({ ...
                {name, 'att', {'char','string'}}
                {A, 'att', 'numeric', 'square'}
                {B, 'att', 'numeric', 'matrix'}
                {c, 'att', 'numeric'}
                {C, 'att', 'numeric', 'matrix'}
                {D, 'att', 'numeric', 'matrix'}
                {k, 'att', 'numeric'}
                {dt, 'att', 'numeric', {'scalar','positive'}}
            });

        end

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

        % check if dimensions fit
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

        if ~isempty(C)
            if ~isscalar(C)
                if size(A,1) ~= size(C,2)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Column dimension of output matrix C must match '...
                        'row dimension of state matrix A.']));
                end
                outputs = size(C,1);
            else % isscalar(C)
                outputs = size(A,1);
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
        
    end

end

function [name,A,B,c,C,D,k,dt,states,inputs,outputs] = ...
    aux_computeProperties(name,A,B,c,C,D,k,dt)
% assign zero vectors/matrices for [] values (from user or by default)
% compute number of states, inputs, and outputs
% difficulty: MATLAB can use 1 instead of eye(n) for multiplication

    % number of states
    states = size(A,1);

    % input matrix
    if ~isempty(A) && isempty(B)
        B = zeros(states,1);
    end

    % number of inputs
    inputs = states;
    if ~isscalar(B)
        inputs = size(B,2);
    end

    % constant offset
    if isempty(c)
        c = zeros(states,1);
    end

    % output matrix
    if isempty(C)
        C = 1;
    end

    % number of outputs
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

end

% ------------------------------ END OF CODE ------------------------------
