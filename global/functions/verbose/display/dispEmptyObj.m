function dispEmptyObj(obj,argname)
% dispEmptyObj - displays text for empty objects
%
% Syntax:
%    dispEmptyObj(obj)
%    dispEmptyObj(obj,argname)
%
% Inputs:
%    obj - object of some class
%    argname - name of obj in workspace
%
% Outputs:
%    -
%
% Example:
%    I = interval();
%    disp(I); % calls dispEmptyObj
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: display

% Authors:       Mark Wetzlinger
% Written:       01-May-2020
% Last update:   02-May-2020
%                05-April-2023 (MW, correct dimensions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp(argname + " =");
fprintf(newline);
n = size(obj);
disp("  " + n(1) + "x" + n(2) + " empty " + class(obj));
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
