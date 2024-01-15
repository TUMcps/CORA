function display(hs)
% display - Displays the properties of a halfspace object (normal vector
%    and distance to the origin) on the command window
%
% Syntax:
%    display(hs)
%
% Inputs:
%    hs - halfspace object
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

% Authors:       Matthias Althoff
% Written:       06-June-2011
% Last update:   02-May-2020 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(hs,'emptySet')
    dispEmptySet(hs,inputname(1));
    return
elseif representsa(hs,'fullspace')
    dispRn(hs,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display normal vector
disp('normal vector: ');
disp(hs.c);

%display distance to origin
disp('distance to origin: ');
disp(hs.d);

% ------------------------------ END OF CODE ------------------------------
