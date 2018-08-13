function display( obj )
% display - Displays the center, generators and constrains of a c-zonotope
%
% Syntax:  
%    display( obj )
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    ---
%
% Example: 
%    Z=contZonotope(rand(2,6));
%    display(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Dmitry Grebenyuk
% Written: 20-December-2017
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

name = [inputname(1), ' = '];
disp(name)

display@zonotope(obj)

%display constraints
disp('A: ');
disp(obj.A);

disp('b: ');
disp(obj.b);

%------------- END OF CODE --------------