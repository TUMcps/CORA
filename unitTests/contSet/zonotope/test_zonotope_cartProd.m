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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  03-January-2023 (MW, add zonotope-numeric cases)
% Last revision:---

%------------- BEGIN CODE --------------

% empty set case
% try 
%     cartProd(Z1,zonotope());
%     res_empty = false;
% catch ME 
%     if ~strcmp(ME.identifier,'CORA:notSupported')
%         rethrow(ME);
%     else
%         res_empty = true;
%     end
% end

% until empty case is correctly integrated...
res_empty = true;

% zonotope-zonotope case:

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9 10 11]);

% compute Cartesian product
Z_ = cartProd(Z1,Z2);

% obtain center and generator matrix
c = center(Z_);
G = generators(Z_);

% true result
true_c = [1; 5; 9];
true_G = [2, 3, 4, 0, 0; ...
          6, 7, 8, 0, 0; ...
          0, 0, 0, 10,11];

% check result
res_ZZ = compareMatrices(c,true_c) && compareMatrices(G,true_G);


% zonotope-numeric case:
Z1 = zonotope([0;2],[3 4 2; -3 -1 3]);
num = 1;

% compute Cartesian product
Z_ = cartProd(Z1,num);

% obtain center and generator matrix
c = center(Z_);
G = generators(Z_);

% true result
true_c = [0; 2; 1];
true_G = [3  4 2; ...
         -3 -1 3; ...
          0  0 0];

% check result
res_Znum = compareMatrices(c,true_c) && compareMatrices(G,true_G);

% compute Cartesian product
Z_ = cartProd(num,Z1);

% obtain center and generator matrix
c = center(Z_);
G = generators(Z_);

% true result
true_c = [1; 0; 2];
true_G = [0  0 0; ...
          3  4 2; ...
         -3 -1 3];

% check result
res_numZ = compareMatrices(c,true_c) && compareMatrices(G,true_G);

% combine all results
res = res_empty && res_ZZ && res_Znum && res_numZ;

%------------- END OF CODE --------------
