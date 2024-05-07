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

% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(fs);
fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
