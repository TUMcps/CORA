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
%    (to console)
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

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               27-Aug-2019
%               01-May-2020 (MW, add empty case)
%               09-June-2020 (MW, restrict number of shown generators)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(Z)
    
    dispEmptyObj(Z,inputname(1));
    
else
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    %display dimension
    display@contSet(Z);
    fprintf(newline);
    
    %display center
    disp('c: ');
    disp(center(Z));

    %display generators
    G = generators(Z);
    maxGens = 10;
    displayGenerators(G,maxGens,'G');
    
end

%------------- END OF CODE --------------