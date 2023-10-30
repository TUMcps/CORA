function display(S)
% display - Displays the properties of a contSet object (dimension) on the
%    command window
%
% Syntax:
%    display(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    ---
%
% Example: 
%    S = contSet();
%    display(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Victor Gassmann
% Written:       02-May-2007
% Last update:   24-March-2022 (VG, replace property read by function call)
%                26-October-2023 (TL, display class properly)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% display dimension
fprintf('%s:\n', class(S))
disp(['- dimension: ', num2str(dim(S))]);

% ------------------------------ END OF CODE ------------------------------
