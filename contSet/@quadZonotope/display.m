function display(Z)
% display - Displays the center and generators of a quadZonotope
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
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      04-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display id, dimension
display(Z.contSet);

%display center
disp('c: ');
disp(Z.c);

%display generators
disp('standard generators: ');
disp(Z.G); 

%display generators
disp('square generators: ');
disp(Z.Gsquare); 

%display generators
disp('quadratic generators: ');
disp(Z.Gquad); 

%display generators
disp('remaining generators: ');
disp(Z.Grest); 

%------------- END OF CODE --------------