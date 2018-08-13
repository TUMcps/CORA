function display(obj)
% display - Displays the center and generators of a zonotope
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - probabilistic zonotope object
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 03-August-2007 
% Last update: 26-February-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%display id, dimension
display(obj.contSet);

%display center
disp('c: ');
disp(obj.Z(:,1));

%display interval generators
disp('interval g_i: ');
disp(obj.Z(:,2:end)); 

%display probabilistic generators
disp('probabilistic g_i: ');
disp(obj.g); 


%display covariance matrix:
disp('covariance matrix: ');
disp(obj.cov);

%------------- END OF CODE --------------