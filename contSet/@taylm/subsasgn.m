function obj = subsasgn(obj, S, value)
% subsasgn - Overloads the opertor that writes elements, e.g. T(1,2)=value,
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    obj = subsasgn(obj, S, value)
%
% Inputs:
%    obj - taylm object 
%    S - contains information of the type and content of element selections
%    value - value to be written
%
% Outputs:
%    obj - taylm object 
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
    obj = builtin('subsasgn', obj, S, value);
    
end

%------------- END OF CODE --------------