function res = assertNarginConstructor(numValidArgs,n_in)
% assertNarginConstructor - asserts if the number of input arguments is
%    matches any entry in a given list of admissible values
%
% Syntax:
%    assertNarginConstructor(numValidArgs,n_in)
%    res = assertNarginConstructor(numValidArgs,n_in)
%
% Inputs:
%    numValidArgs - ordered list of admissible number of input arguments
%    n_in - number of input arguments
%
% Outputs:
%    res - true if assertion is successful
%
% Example:
%    assertNarginConstructor([0,1,4,5],1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: assert, assertLoop, assertThrowsAs

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different between finite admissible list and infinite
if numValidArgs(end) == Inf
    % if list of admissible numbers of input arguments ends in Inf, any
    % value between the penultimate value and Inf is also allowed
    if ~any(numValidArgs(1:end-1) == n_in) && n_in < numValidArgs(end-1)
        throwAsCaller(CORAerror('CORA:numInputArgsConstructor',numValidArgs));
    end
else
    % no Inf -> must be any of the given values
    if ~any(numValidArgs == n_in)
        throwAsCaller(CORAerror('CORA:numInputArgsConstructor',numValidArgs));
    end
end

if nargout == 1
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
