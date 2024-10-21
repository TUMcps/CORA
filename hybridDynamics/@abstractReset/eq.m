function res = eq(reset1,reset2)
% eq - checks if two reset functions have equal pre-/post-state and input
%    dimensions
%
% Syntax:
%    res = reset1 == reset2
%    res = eq(reset1,reset2)
%
% Inputs:
%    reset1 - abstractReset object
%    reset2 - abstractReset object
%
% Outputs:
%    res - true/false
%
% Example: 
%    reset1 = abstractReset(2,1,2);
%    reset2 = abstractReset(2,1,3);
%    reset1 == reset1;
%    reset1 == reset2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reset/isequal

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% redirect to isequal
res = isequal(reset1,reset2);

% ------------------------------ END OF CODE ------------------------------
