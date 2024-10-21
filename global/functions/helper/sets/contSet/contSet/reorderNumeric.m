function [S1,S2] = reorderNumeric(S1,S2)
% reorderNumeric - reorders the input arguments so that the first one is a
%    class object and the second one is numeric; this occurs, e.g., if a
%    contSet function is called with numeric type as a first input argument
%    and the respective contSet class as a second input argument
%
% Syntax:
%    [S1,S2] = reorderNumeric(S1,S2)
%
% Inputs:
%    S1 - class object or numeric
%    S2 - class object or numeric
%
% Outputs:
%    S1 - class object
%    S2 - numeric
%
% Example:
%    S1 = 4;
%    S2 = interval(-2,3);
%    [S1,S2] = reorderNumeric(S1,S2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% classic temporary-swap
if isnumeric(S1)
    S1_copy = S1;
    S1 = S2;
    S2 = S1_copy;
end

% ------------------------------ END OF CODE ------------------------------
