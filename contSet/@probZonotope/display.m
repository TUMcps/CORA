function display(probZ)
% display - Displays the properties of a probZonotope object (center,
%    interval generators, probabilistic generators, covariance matrix) on
%    the command window
%
% Syntax:
%    display(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
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

% Authors:       Matthias Althoff
% Written:       03-August-2007 
% Last update:   26-February-2008
%                09-June-2020 (MW, update formatting of output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(probZ);
fprintf(newline);

%display center
disp('center: ');
disp(center(probZ));

maxGens = 10;

%display interval generators
displayGenerators(probZ.Z(:,2:end),maxGens,'interval generators');

%display probabilistic generators
displayGenerators(probZ.g,maxGens,'probabilistic generators');

%display covariance matrix:
disp('covariance matrix: ');
disp(probZ.cov);

% ------------------------------ END OF CODE ------------------------------
