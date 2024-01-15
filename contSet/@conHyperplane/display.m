function display(hyp)
% display - Displays the properties of a conHyperplane object (halfspace,
%    constraint system) on the command window
%
% Syntax:
%    display(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    ---
%
% Example:
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       10-August-2011
% Last update:   02-May-2020 (MW, added empty case)
%                18-June-2022 (MW, empty constraint system)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa_(hyp,'emptySet',0)
    dispEmptySet(hyp,inputname(1));
    return
elseif representsa_(hyp,'fullspace',0)
    dispRn(hyp,inputname(1));
    return
end

try
    inputname;
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
catch
    % nothing here
end

%display hyperplane
disp('normal vector:');
disp(hyp.a);
disp('distance to origin:');
disp(hyp.b);

%display constraint system
if representsa_(hyp,'hyperplane',eps)
    disp('constraint system (Cx <= d): (none)');
else
    disp('constraint system (Cx <= d):');

    disp('C:');
    disp(hyp.C);
    disp('d:');
    disp(hyp.d);
end

% ------------------------------ END OF CODE ------------------------------
