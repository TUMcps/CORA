function display(Z)
% display - Displays the center and generators of a zonotope
%
% Syntax:  
%    display(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    ---
%
% Example: 
%    Z=zonotope(rand(2,6));
%    display(Z);
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

%display id, dimension
display(Z.contSet);

%display center
disp('c: ');
disp(Z.Z(:,1));

%display generators
disp('g_i: ');
disp(Z.Z(:,2:end)); 

%------------- END OF CODE --------------