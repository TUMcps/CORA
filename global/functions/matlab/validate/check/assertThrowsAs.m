function res = assertThrowsAs(fun,id,varargin)
% assertThrowsAs - asserts if the execution of a function handle throws a 
%    specific error when called with the arguments in varargin
%
% Syntax:
%    assertThrowsAs(fun,id,varargin)
%    res = assertThrowsAs(fun,id,varargin)
%
% Inputs:
%    fun - function handle
%    id - MException.identifier
%    varargin - variable number of input arguments to function handle
%
% Outputs:
%    res - true if assertion is successful
%
% Example:
%    assertThrowsAs(@capsule,'CORA:wrongInputInConstructor',1,[2;1],0.5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: assert, assertLoop

% Authors:       Mark Wetzlinger
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    fun(varargin{:});
    assert(false,sprintf("Function call should have thrown the error '%s' but threw nothing.", id));
catch ME
    if ~strcmp(ME.identifier,id)
        rethrow(ME);
    end
    res = true;
end

% no output if not desired
if nargout == 0
    clear res
end

% ------------------------------ END OF CODE ------------------------------
