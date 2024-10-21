classdef abstractReset
% abstractReset - abstract superclass for reset functions
%
% Syntax:
%    reset = abstractReset(preStateDim,inputDim,postStateDim)
%
% Inputs:
%    preStateDim - dimension of state before reset
%    inputDim - dimension of input before reset
%    postStateDim - dimension of state after reset
%
% Outputs:
%    reset - generated abstractReset object
%
% Example:
%    reset = abstractReset(2,1,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    
    preStateDim;        % state dimension before reset
    inputDim;           % input dimension
    postStateDim;       % state dimension after reset

end

methods
    
    % class constructor
    function reset = abstractReset(varargin)

        % 0. check number of input arguments
        assertNarginConstructor([1,3],nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'abstractReset')
            reset = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [preStateDim,inputDim,postStateDim] = setDefaultValues({0,0,0},varargin);

        % 3. check correctness of input arguments
        if CHECKS_ENABLED && nargin > 0
            inputArgsCheck({ ...
                {preStateDim, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {inputDim, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {postStateDim, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}};
            })
        end

        % 4. assign properties
        reset.preStateDim = preStateDim;
        reset.inputDim = inputDim;
        reset.postStateDim = postStateDim;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
