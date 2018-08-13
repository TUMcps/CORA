function obj = removeZeroConstraints(obj)
% removeZeroConstraints - remove the trivial constraint [0 0 .. 0]*ksi = 0
%
% Syntax:  
%    obj = removeZeroConstraints(obj)
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    obj - c-zonotope object without trivial constraints
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      29-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Remove the trivial constraint [0 0 .. 0]*ksi = 0

    % select all trivial constraints
    temp = all(obj.A==0,2);
    ind = find(temp == 1);
    obj.A(ind,:) = [];
    
    % check if all eliminated constraints are satisfiable
    if any(obj.b(ind))
       error('Unsatisfiable constraint [0 0 ...0]*x = c !'); 
    end

    obj.b(ind) = [];
    
end

%------------- END OF CODE --------------