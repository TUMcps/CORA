function dispEmptySet(S,argname)
% dispEmptySet - displays text for empty objects
%
% Syntax:
%    dispEmptySet(S)
%    dispEmptySet(S,argname)
%
% Inputs:
%    S - contSet object
%    argname - name of obj in workspace
%
% Outputs:
%    -
%
% Example:
%    I = interval.empty(2)
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
disp("  " + dim(S) + "-dimensional empty set (represented as " + class(S) + ")");
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
