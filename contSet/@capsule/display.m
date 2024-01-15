function display(C)
% display - Displays the properties of a capsule object (center, generator,
%    radius) on the command window
%
% Syntax:
%    display(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    ---
%
% Example: 
%    C = capsule([1;1],[1;0],0.5);
%    display(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2019
% Last update:   02-May-2020 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(C,'emptySet')
    dispEmptySet(C,inputname(1));
    return
elseif representsa(C,'fullspace')
    dispRn(C,inputname(1));
    return
end

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display dimension
display@contSet(C);
fprintf(newline);

% display center
disp('center: ');
disp(C.c);

% display generator
disp('generator: ');
disp(C.g); 

% display radius
disp('radius: ');
disp(C.r);

% ------------------------------ END OF CODE ------------------------------
