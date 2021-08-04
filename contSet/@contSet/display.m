function display(obj)
% display - Displays the id and dimension of a continuous set
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - contSet object
%
% Outputs:
%    ---
%
% Example: 
%    S=contSet(2);
%    display(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       02-May-2007
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% display dimension
disp(['dimension: ', num2str(obj.dimension)]);

%------------- END OF CODE --------------