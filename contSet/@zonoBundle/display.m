function display(obj)
% display - Displays a zonotope bundle
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    Z - zonotope bundle
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
% Written:      09-November-2010
% Last update:  02-May-2020 (MW, add empty case)
%               09-June-2020 (MW, remove dependency from zonotope/display)
% Last revision:---

%------------- BEGIN CODE --------------

%display each zonotope
if obj.parallelSets == 0
    dispEmptyObj(obj,inputname(1));
    
else
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    % display id, dimension
    display@contSet(obj);
    fprintf(newline);
    
    maxGens = 10;
    for i=1:obj.parallelSets
        
        disp(['zonotope ',num2str(i),':',newline]);
        
        %display center
        disp('c: ');
        disp(center(obj.Z{i}));

        %display generators
        G = generators(obj.Z{i});
        displayGenerators(G,maxGens,'G');
        
    end
    
end

%------------- END OF CODE --------------