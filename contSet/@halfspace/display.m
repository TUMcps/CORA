function display(obj)
% display - Displays the normal vector and distance to the origin of a
% halfspace
%
% Syntax:  
%    display(h)
%
% Inputs:
%    h - halfspace object
%
% Outputs:
%    ---
%
% Example: 
%    h = halfspace([1;2],3);
%    display(h)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-June-2011
% Last update:  02-May-2020 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(obj)
    
    dispEmptyObj(obj,inputname(1));
    
else

    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);

    %display normal vector
    disp('normal vector: ');
    disp(obj.c);

    %display distance to origin
    disp('distance to origin: ');
    disp(obj.d);

end

%------------- END OF CODE --------------