function display(V)
% display - Displays the matrix of vertices as column vectors
%
% Syntax:  
%    display(V)
%
% Inputs:
%    V - vertices object
%
% Outputs:
%    ---
%
% Example: 
%    V=vertices(rand(2,6));
%    display(V);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 14-September-2006 
% Last update: 22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

disp('V: ');
disp(V.V);

%------------- END OF CODE --------------