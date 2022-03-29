function res = levelSet(obj)
% levelSet - converts mptPolytope object to equivalent levelSet object
%
% Syntax:  
%    res = levelSet(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - levelSet object representing the same set as obj
%
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author:       Maximilian Perschl
% Written:      23-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

vars = sym('x',[dim(obj), 1]);

A = obj.P.A;
b = obj.P.b;

if ~isempty(obj.P.Ae)
    error("Not implemented yet!");
end

equations = A*vars - b;
compOps = repmat({'<='},size(A,1),1);

res = levelSet(equations,vars,compOps);

  
end

%------------- END OF CODE --------------
