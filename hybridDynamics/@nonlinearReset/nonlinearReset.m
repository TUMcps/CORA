classdef nonlinearReset < abstractReset
% nonlinearReset - constructor for nonlinear reset functions
%
% Description:
%    This class represents nonlinear reset functions
%        x_ = f(x,u)
%    where x and u are the state and input before transition, and x_ is the
%    state after transition
%
% Syntax:
%    nonlinReset = nonlinearReset()
%    nonlinReset = nonlinearReset(f)
%    nonlinReset = nonlinearReset(f,preStateDim,inputDim,postStateDim)
%
% Inputs:
%    f - function handle
%    preStateDim - length of vector x
%    inputDim - length of vector u
%    postStateDim - length of vector x_
%
% Outputs:
%    nonlinReset - generated nonlinearReset object
%
% Example:
%    f = @(x,u) [x(1) + 2*u(1); x(2)];
%    nonlinReset = nonlinearReset(f);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reset, linearReset

% Authors:       Mark Wetzlinger
% Written:       07-September-2024
% Last update:   15-October-2024 (MW, add dimensions to input arguments)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    f;                  % function handle to mapping matrix

    % --- computed by nonlinearReset/derivatives.m ---
    tensorOrder;        % tensor order for evaluation
    J;                  % function handle to Jacobian matrix
    H;                  % function handle to Hessian matrix
    T;                  % function handle to third-order tensor
end

methods
    
    % class constructor
    function nonlinReset = nonlinearReset(varargin)

        % 0. check number of input arguments
        assertNarginConstructor([0,1,4],nargin);

        % 1. copy constructor: not allowed due to obj@abstractReset below...
        % if nargin == 1 && isa(varargin{1},'nonlinearReset')
        %     nonlinReset = varargin{1}; return
        % end
            
        % 2. parse input arguments: varargin -> vars
        [f,preStateDim,inputDim,postStateDim] = aux_parseInputArgs(varargin{:});
        
        % 3. check correctness of input arguments
        aux_checkInputArgs(f,preStateDim,inputDim,postStateDim,nargin);

        % 4. compute number of states, inputs, and outputs
        [f,preStateDim,inputDim,postStateDim] = ...
            aux_computeProperties(f,preStateDim,inputDim,postStateDim);
        
        % 5. instantiate parent class, assign properties
        nonlinReset@abstractReset(preStateDim,inputDim,postStateDim); 
        nonlinReset.f = f;

    end

end

end


% Auxiliary functions -----------------------------------------------------

function [f,preStateDim,inputDim,postStateDim] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % init properties
    f = @(x,u) []; preStateDim = []; inputDim = []; postStateDim = [];

    % no input arguments
    if nargin == 0
        return
    end

    % only function handle given
    if nargin == 1
        f = varargin{1};
        return
    end

    % all dimensions must be given
    narginchk(4,4);

    % set defaults
    [f,preStateDim,inputDim,postStateDim] = ...
        setDefaultValues({f,preStateDim,inputDim,postStateDim},varargin);
    
end

function aux_checkInputArgs(f,preStateDim,inputDim,postStateDim,n_in)
% ensure that nonlinear reset function x_ = f(x,u) is properly defined

    if CHECKS_ENABLED && n_in > 0
        
        inputArgsCheck({ ...
            {f, 'att', 'function_handle'} ...
            {preStateDim, 'att', 'numeric'} ...
            {inputDim, 'att', 'numeric'} ...
            {postStateDim, 'att', 'numeric'}
        });

        % if preStateDim and inputDim provided, check if they are plausible
        % by inserting vector of given size into the function handle
        if ~isempty(preStateDim) && ~isempty(inputDim)
            try
                f_out = f(zeros(preStateDim,1),zeros(inputDim,1));
            catch ME
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Sizes for pre-state and input do not match the given function handle.'));
            end
            % check if output has appropriate size (must be vector)
            if ~isempty(postStateDim)
                if ~iscolumn(f_out)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Function handle must return a column vector.'));
                elseif size(f_out,1) ~= postStateDim
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Size for post-state does not match the given function handle.'));
                end
            end
        end        
    end

end

function [f,preStateDim,inputDim,postStateDim] = aux_computeProperties(f,preStateDim,inputDim,postStateDim)

    % as long as aux_parseInputArgs enforces either 1 or 4 input arguments,
    % we must only check 'preStateDim' and can then compute all values
    if isempty(preStateDim)
        [count,postStateDim] = inputArgsLength(f,2);
        preStateDim = count(1);
        % for ease of code, we always have at least one input dimension
        inputDim = count(2);
    end

    % input dimension must also be >= 1
    inputDim = max(inputDim,1);

end

% ------------------------------ END OF CODE ------------------------------
