function [qZ]=reduceGirard(qZ,order)
% reduceGirard - Reduce remaining generators of a quadratic zonotope so 
% that its order stays below a specified limit 
%
% Syntax:  
%    [qZred]=reduceGirard(qZ,order)
%
% Inputs:
%    qZ - quadZonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    qZ - reduced quadZonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      05-September-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%dim
dim = length(qZ.c);

%create zonotope of the remaining generators
Z = zonotope([zeros(dim,1),qZ.Grest]);

%initialize Z_red
Zred = reduce(Z,'girard',order);

%obtain reduced generators
Zmat = get(Zred,'Z');
qZ.Grest = Zmat(:,2:end);

%------------- END OF CODE --------------