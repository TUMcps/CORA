function display(E)
% display - Displays the properties of an ellipsoid object (center, shape
%    matrix, dimension, degeneracy) on the command window
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
%    E = ellipsoid([1,0;0,1]);
%    display(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   02-May-2020 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'}});

% special cases
if representsa(E,'emptySet')
    dispEmptySet(E,inputname(1));
    return
elseif representsa(E,'fullspace')
    dispRn(E,inputname(1));
    return
end


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
disp(dim(E)); 

%display whether degenerate or not
disp('degenerate: ');
disp(~isFullDim(E));

% ------------------------------ END OF CODE ------------------------------
