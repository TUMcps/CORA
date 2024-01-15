function dispRn(S,argname)
% dispRn - displays text for fullspace objects
%
% Syntax:
%    dispRn(S)
%    dispRn(S,argname)
%
% Inputs:
%    S - contSet object
%    argname - name of obj in workspace
%
% Outputs:
%    -
%
% Example:
%    I = interval.Inf(2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: display

% Authors:       Mark Wetzlinger
% Written:       10-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp(argname + " =");
fprintf(newline);

% display
disp("  R^" + dim(S) + " (represented as " + class(S) + ")");
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
