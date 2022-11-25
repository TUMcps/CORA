function display(obj)
% display - Displays the center, generators and constraints of a conZonotope
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    ---
%
% Example: 
%    Z = conZonotope(rand(2,6));
%    display(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Dmitry Grebenyuk, Mark Wetzlinger
% Written:       20-December-2017
% Last update:   01-May-2020 (MW, added empty case)
%                09-June-2020 (MW, restrict number of shown generators)
% Last revision: ---

%------------- BEGIN CODE --------------

if isempty(obj)
    
    dispEmptyObj(obj,inputname(1));

else
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);

    % display center and generators
    disp('c: ');
    disp(obj.Z(:,1));

    G = obj.Z(:,2:end);
    maxGens = 10;
    displayGenerators(G,maxGens,'G');
    
    %display constraint system
    if isempty(obj.A) && isempty(obj.b)
        disp('Constraint system (Cx <= d): no constraints.')
        fprintf(newline);
        
    else
        disp('Constraint system (Cx <= d):')

        displayGenerators(obj.A,maxGens,'A');

        disp('b: ');
        disp(obj.b);
    end
    
end


end
%------------- END OF CODE --------------