function display(E)
% display - Displays the center and generators of an Ellipsoid
%
% Syntax:  
%    display(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    ---
%
% Example: 
%    E=ellipsoid([1,0;0,1]);
%    display(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  02-May-2020 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(E)
    
    dispEmptyObj(E,inputname(1));
    
else

    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    %display center
    disp('q: ');
    disp(E.q);

    %display shape matrix
    disp('Q: ');
    disp(E.Q); 

    %display actual dimension
    disp('dimension: ');
    disp(E.dim); 

    %display whether degenerate or not
    disp('degenerate: ');
    disp(E.isdegenerate);

end

%------------- END OF CODE --------------