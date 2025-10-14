classdef linearSysDT < contDynamics
% linearSysDT - object constructor for linear discrete-time systems
%    
% Description:
%    Generates a discrete-time linear system object according to the 
%    following first-order difference equations:
%       x(k+1) = A x(k) + B u(k) + c + E w(k)
%       y(k)   = C x(k) + D u(k) + k + F v(k)
%
% Syntax:
%    linsysDT = linearSysDT(A,B,dt)
%    linsysDT = linearSysDT(A,B,c,dt)
%    linsysDT = linearSysDT(A,B,c,C,dt)
%    linsysDT = linearSysDT(A,B,c,C,D,dt)
%    linsysDT = linearSysDT(A,B,c,C,D,k,dt)
%    linsysDT = linearSysDT(A,B,c,C,D,k,E,dt)
%    linsysDT = linearSysDT(A,B,c,C,D,k,E,F,dt)
%    linsysDT = linearSysDT(name,A,B,dt)
%    linsysDT = linearSysDT(name,A,B,c,dt)
%    linsysDT = linearSysDT(name,A,B,c,C,dt)
%    linsysDT = linearSysDT(name,A,B,c,C,D,dt)
%    linsysDT = linearSysDT(name,A,B,c,C,D,k,dt)
%    linsysDT = linearSysDT(name,A,B,c,C,D,k,E,dt)
%    linsysDT = linearSysDT(name,A,B,c,C,D,k,E,F,dt)
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
%    dt - sampling time
%
% Outputs:
%    linsysDT - generated linearSysDT object
%
% Example:
%    A = [-0.4 0.6; 0.6 -0.4];
%    B = [0; 1];
%    c = [0;0];
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
%                30-August-2024 (MW, integrate E and F matrices)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A   % system matrix: n x n
    B   % input matrix: n x m
    c   % constant input: n x 1
    C   % output matrix: q x n
    D   % feedthrough matrix: q x m
    k   % output offset: q x 1
    E   % disturbance matrix: n x r
    F   % output disturbance matrix: q x s
    dt  % sampling time
end

methods
    
    % class constructor
    function linsysDT = linearSysDT(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:10,nargin);

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'linearSysDT')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,A,B,c,C,D,k,E,F,dt] = aux_parseInputArgs(varargin{:});

        if iscell(A) && iscell(B) && iscell(C) && iscell(D)
            % --> linear time-varying system

            if ~isempty(c) || ~isempty(k) || ~isempty(E) || ~isempty(F)
                throw(CORAerror('CORA:notSupported',...
                    "Support for LTV systems with nonempty c, k, E or F is not implemented yet"));
            end
            for i=1:length(A)
                % 3. check correctness of input arguments
                aux_checkInputArgs(name,A{i},B{i},c,C{i},D{i},k,E,F,dt,nargin);

                % 4. compute number of states, inputs, and outputs
                [name,A{i},B{i},c,C{i},D{i},k,E,F,dt,states,inputs,outputs,dists,noises] = ...
                    aux_computeProperties(name,A{i},B{i},c,C{i},D{i},k,E,F,dt);
            end
        else
            % 3. check correctness of input arguments
            aux_checkInputArgs(name,A,B,c,C,D,k,E,F,dt,nargin);

            % 4. compute number of states, inputs, and outputs
            [name,A,B,c,C,D,k,E,F,dt,states,inputs,outputs,dists,noises] = ...
                aux_computeProperties(name,A,B,c,C,D,k,E,F,dt);
        end
        
        % 5. instantiate parent class, assign properties
        linsysDT@contDynamics(name,states,inputs,outputs,dists,noises); 
        linsysDT.A = A; linsysDT.B = B; linsysDT.c = c;
        linsysDT.C = C; linsysDT.D = D; linsysDT.k = k;
        linsysDT.E = E; linsysDT.F = F;
        linsysDT.dt = dt;
    end
    
    % invoke function observe so that the superclass can access private
    % functions
    function [R, tcomp] = observe(syslinDT,params,options)
        [R, tcomp] = observe@contDynamics(syslinDT,params,options);
    end

    % update system dynamics for the new augmented input [u; w] where w is
    % the process noise acting on all states 
    function syslinDT = augment_u_with_w(syslinDT)
        if ~isempty(syslinDT.E)
            E = syslinDT.E;
        else
            E = eye(syslinDT.nrOfDims);
        end
        syslinDT = linearSysDT(syslinDT.name,syslinDT.A, [syslinDT.B E], [], ...
            syslinDT.C, [syslinDT.D, zeros(syslinDT.nrOfOutputs, size(E,2))], [],[],syslinDT.F, syslinDT.dt);
    end

    % update system dynamics for the new augmented input [u; v] where v is
    % the measurement noise acting on all outputs 
    function syslinDT = augment_u_with_v(syslinDT)
        if ~isempty(syslinDT.E)
            F = syslinDT.F;
        else
            F = eye(syslinDT.nrOfOutputs);
        end
        syslinDT = linearSysDT(syslinDT.A, [syslinDT.B zeros(syslinDT.nrOfDims, size(F,2))], [], ...
            syslinDT.C, [syslinDT.D, F], [], syslinDT.E, [], syslinDT.dt);
    end
end

methods (Static = true)
    linsys = identify(varargin)       % identifies linear system from data
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,c,C,D,k,E,F,dt] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % default name
    def_name = 'linearSysDT';
    A = []; B = []; c = []; C = []; D = []; k = []; E = []; F = []; dt = 0;

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
        [name,A,B,c,C,D,k,E,F] = setDefaultValues({def_name,A,B,c,C,D,k,E,F},varargin);

    else
        % set default name
        name = def_name;
        [A,B,c,C,D,k,E,F] = setDefaultValues({A,B,c,C,D,k,E,F},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,c,C,D,k,E,F,dt,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % ensure that values have correct data type
        nameCheck = {};
        if ~strcmp(name,'linearSysDT')
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
            {dt, 'att', 'numeric', {'scalar','positive'}}
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

        % check input matrix and throughput matrix
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

        % check output matrix and throughput matrix
        outputs = size(A,1); % assume outputs = states
        if ~isempty(C)
            if ~isscalar(C)
                if size(A,1) ~= size(C,2)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Column dimension of output matrix C must match '...
                        'row dimension of state matrix A.']));
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

        % check if output disturbance matrix has the correct size
        if ~isempty(F)
            if ~isscalar(F) && outputs ~= size(F,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Column dimension of noise matrix F must match '...
                    'row/column dimension of state matrix C.']));
            end
        end
        
    end

end

function [name,A,B,c,C,D,k,E,F,dt,states,inputs,outputs,dists,noises] = ...
    aux_computeProperties(name,A,B,c,C,D,k,E,F,dt)
% assign zero vectors/matrices for [] values (from user or by default)
% compute number of states, inputs,  outputs, disturbances, and noises
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

    % disturbance matrix and number of disturbances
    if isempty(E)
        E = zeros(states,1);
    end
    dists = states;
    if ~isscalar(E)
        dists = size(E,2);
    end

    % noise matrix and number of noises
    if isempty(F)
        F = zeros(outputs,1);
    end
    noises = outputs;
    if ~isscalar(F)
        noises = size(F,2);
    end

end

% ------------------------------ END OF CODE ------------------------------
