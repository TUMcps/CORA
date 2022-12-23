function res = test_zonotope_cartProd
% test_zonotope_cartProd - unit test function of cartesian product
%
% Syntax:  
%    res = test_zonotope_cartProd
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9 10 11]);

% empty set
% try 
%     cartProd(Z1,zonotope());
%     res_e = false;
% catch ME 
%     if ~strcmp(ME.identifier,'CORA:notSupported')
%         rethrow(ME);
%     else
%         res_e = true;
%     end
% end

% compute Cartesian product
Z3 = cartProd(Z1,Z2);

% obtain zonotope matrix
c = center(Z3);
G = generators(Z3);

% true result
true_c = [1; 5; 9];
true_G = [2, 3, 4, 0, 0; ...
          6, 7, 8, 0, 0; ...
          0, 0, 0, 10,11];

% check result
res = compareMatrices(c,true_c) && compareMatrices(G,true_G); % && res_e;

%------------- END OF CODE --------------
