classdef linProbSys < contDynamics
% linProbSys - class constructor for linear probabilistic systems
%
% Description:
%    Generates a linear stochastic differential system, also known as the
%    multivariate Ornstein-Uhlenbeck process:
%       x'(t) = A x(t) + B u(t) + C xi(t),
%    where xi(t) is white Gaussian noise.
%
% Syntax:
%    obj = linProbSys(A,B)
%    obj = linProbSys(A,B,C)
%    obj = linProbSys(name,A,B)
%    obj = linProbSys(name,A,B,C)
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    C - noise matrix
%
% Outputs:
%    obj - generated linProbSys object
%
% Example:
%    A = [-1 -4; 4 -1];
%    B = eye(2);
%    C = 0.7*eye(2);
%    sys = linProbSys(A,B,C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       06-October-2007 
% Last update:   26-February-2008
%                05-August-2016 (changed to new OO format)
%                19-June-2022 (MW, update syntax)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A;      % state matrix
    B;      % input matrix
    C;      % noise matrix

    % internally-set properties
    taylor = [];
end

methods
    % class constructor
    function obj = linProbSys(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'linProbSys')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,A,B,C] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,A,B,C,nargin);

        % 4. compute number of states and inputs
        [name,A,B,C,states,inputs] = aux_computeProperties(name,A,B,C);
        
        % 5. instantiate parent class, assign properties
        obj@contDynamics(name,states,inputs,0); 
        obj.A = A; obj.B = B; obj.C = C;

    end
end
end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,C] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    end

    % default name
    def_name = 'linProbSys';
    % init properties
    A = []; B = []; C = [];

    % no input arguments
    if nargin == 0
        name = def_name;
        return
    end

    % parse depending on whether first input argument is the name
    if ischar(varargin{1})
        % first input argument: name
        [name,A,B,C] = setDefaultValues({def_name,A,B,C},varargin);

    else
        % set default name
        name = def_name;
        [A,B,C] = setDefaultValues({A,B,C},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,C,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % ensure that values have correct data type
        if strcmp(name,'linProbSys')
            % default name (unless explicitly chosen by user, we have A as
            % first input argument)

            inputArgsCheck({ ...
                {A, 'att', 'numeric', 'square'}
                {B, 'att', 'numeric', 'matrix'}
                {C, 'att', 'numeric', 'square'}
            });
            
        else

            inputArgsCheck({ ...
                {name, 'att', {'char','string'}}
                {A, 'att', 'numeric', 'square'}
                {B, 'att', 'numeric', 'matrix'}
                {C, 'att', 'numeric', 'square'}
            });

        end

        % check if dimensions fit
        if ~isempty(B) && ~isscalar(B) && size(A,1) ~= size(B,1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
            ['Column dimension of input matrix B must match '...
            'row/column dimension of state matrix A.']));
        end
        if ~isempty(C) && ~isscalar(C) && size(A,1) ~= size(C,2)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Row/column dimension of noise matrix C must match '...
                'row/column dimension of state matrix A.']));
        end
        
    end

end

function [name,A,B,C,states,inputs] = aux_computeProperties(name,A,B,C)
% assign zero vectors/matrices for [] values (from user or by default)
% compute number of states and inputs
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

    % noise matrix
    if isempty(C)
        C = 1;
    end

end

% ------------------------------ END OF CODE ------------------------------
