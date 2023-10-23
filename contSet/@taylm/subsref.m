function res = subsref(tay,S)
% subsref - Overloads the operator that selects elements, e.g. T(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:
%    res = subsref(tay,S)
%
% Inputs:
%    tay - a taylm object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    res - element or elemets of the taylm matrix
%
% Example:
%    tay = taylm(interval([1;3],[2;4]),4,{'x';'y'});
%    tay(1)
%    tay(2)
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Authors:       Dmitry Grebenyuk
% Written:       20-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call built-in function
res = builtin('subsref', tay, S);

% ------------------------------ END OF CODE ------------------------------
