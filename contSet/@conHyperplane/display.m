function display(obj)
% display - Displays the hyperplane and the constraint system
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - conHyperplane object
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

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  02-May-2020 (MW, added empty case)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(obj)
    
    dispEmptyObj(obj,inputname(1));
    
else

    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    %display hyperplane
    disp('normal vector:');
    disp(obj.h.c);
    disp('distance to origin:');
    disp(obj.h.d);

    %display constraint system
    disp('Constraint system (Cx <= d):')

    disp('C:');
    disp(obj.C);
    disp('d:');
    disp(obj.d);

end

%------------- END OF CODE --------------