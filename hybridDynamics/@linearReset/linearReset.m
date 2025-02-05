classdef linearReset < abstractReset
% linearReset - constructor for linear reset functions
%
% Description:
%    This class represents linear reset functions
%        x_ = Ax + Bu + c
%    where x and u are the state and input before transition, and x_ is the
%    state after transition
%
% Syntax:
%    linReset = linearReset()
%    linReset = linearReset(A)
%    linReset = linearReset(A,B)
%    linReset = linearReset(A,B,c)
%
% Inputs:
%    A - state mapping matrix
%    B - input mapping matrix
%    c - constant translation
%
% Outputs:
%    linReset - generated linearReset object
%
% Example:
%    A = [1 2; -1 1; 0 1];
%    B = [1; -1; 0];
%    c = [-1; 0; 0];
%
%    linReset = linearReset(A,B,c);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reset, nonlinearReset

% Authors:       Mark Wetzlinger
% Written:       07-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A;  % state mapping matrix
    B;  % input mapping matrix
    c;  % constant offset
end

methods
    
    % class constructor
    function linReset = linearReset(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:3,nargin);

        % 1. copy constructor: not allowed due to obj@abstractReset below
%         if nargin == 1 && isa(varargin{1},'linearReset')
%             linReset = varargin{1}; return
%         end
            
        % 2. parse input arguments: varargin -> vars
        [A,B,c] = aux_parseInputArgs(varargin{:});
        
        % 3. check correctness of input arguments
        aux_checkInputArgs(A,B,c,nargin);

        % 4. compute number of states, inputs, and outputs
        [A,B,c,preStateDim,inputDim,postStateDim] = ...
            aux_computeProperties(A,B,c);
        
        % 5. instantiate parent class, assign properties
        linReset@abstractReset(preStateDim,inputDim,postStateDim); 
        linReset.A = A; linReset.B = B; linReset.c = c;

    end

end

methods (Static = true)
    linReset = eye(varargin); % identity reset function of dimension n
end

end


% Auxiliary functions -----------------------------------------------------

function [A,B,c] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % init properties
    A = []; B = []; c = [];

    % no input arguments
    if nargin == 0
        return
    end

    % set defaults
    [A,B,c] = setDefaultValues({A,B,c},varargin);
    
end

function aux_checkInputArgs(A,B,c,n_in)
% ensure that linear reset function x_ = Ax + Bu + c is properly defined

    if CHECKS_ENABLED && n_in > 0
        
        inputArgsCheck({ ...
            {A, 'att', 'numeric', 'matrix'}
            {B, 'att', 'numeric'}
            {c, 'att', 'numeric'}
        });

        % check size
        if ~isempty(c) && ~isvector(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Offset c must be a vector.'));
        end
        if ~isempty(c) && ~isempty(A) && size(A,1) ~= length(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Length of offset c must match row dimension of state mapping matrix A.'));
        end
        if ~isempty(B) && ~isempty(A) && size(A,1) ~= size(B,1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Row dimension of input mapping matrix B must match row dimension of state mapping matrix A.'));
        end
        
    end

end

function [A,B,c,preStateDim,inputDim,postStateDim] = aux_computeProperties(A,B,c)
    % compute properties

    if isempty(A)
        if ~isempty(B)
            A = zeros(size(B,1),0);
        elseif ~isempty(c)
            A = zeros(size(c,1),0);
        end
    end

    [postStateDim, preStateDim] = size(A);
    
    inputDim = size(B,2);
    if isempty(B)
        B = zeros(postStateDim,1);
        inputDim = 1;
    end

    if isempty(c)
        c = zeros(postStateDim,1);
    end

end

% ------------------------------ END OF CODE ------------------------------
