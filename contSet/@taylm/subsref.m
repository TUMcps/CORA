function res = subsref(obj, S)
% subsref - Overloads the opertor that selects elements, e.g. T(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    res = subsref(obj, S)
%
% Inputs:
%    obj - a taylm object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    res - element or elemets of the taylm matrix
%
% Example: 
%    t(1,2).monomials
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      20-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % call built-in function
    res = builtin('subsref', obj, S);
    
end
%------------- END OF CODE --------------
