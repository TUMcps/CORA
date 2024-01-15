function display(fs)
% display - displays the properties of a fullspace object (dimension) on
%    the command window
%
% Syntax:
%    display(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    -
%
% Example: 
%    fs = fullspace(2);
%    display(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display R^n
disp("  R^" + fs.dimension);
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
