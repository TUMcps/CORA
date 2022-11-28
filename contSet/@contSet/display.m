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
%    S = contSet(2);
%    display(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Victor Gassmann
% Written:       02-May-2007
% Last update:   24-March-2022 (VG: replace property read by function call)
% Last revision: ---

%------------- BEGIN CODE --------------

% display dimension
disp(['dimension: ', num2str(dim(S))]);

%------------- END OF CODE --------------