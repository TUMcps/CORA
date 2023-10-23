function display(O)
% display - displays the properties of an emptySet object (dimension) on
%    the command window
%
% Syntax:
%    display(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    -
%
% Example: 
%    O = emptySet(2);
%    display(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display 
disp("  " + O.dimension + "-dimensional empty set");
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
