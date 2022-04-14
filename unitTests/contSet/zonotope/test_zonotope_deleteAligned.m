function res = test_zonotope_deleteAligned
% test_zonotope_deleteAligned - unit test function of deleteAligned
%
% Syntax:  
%    res = test_zonotope_deleteAligned
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

% Author:       Mark Wetzlinger
% Written:      26-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create a zonotope
Z_cent = zeros(2,1);
% aligned generators differ by scaling, sign
Z_gens = [4 2 2 3 1 -4;
          2 3 1 0 2 -2];
Z = zonotope([Z_cent, Z_gens]);

% delete aligned generators
Z_del = deleteAligned(Z);

% true matrix
Z_gens_true = [10 2 3 1;
               5  3 0 2];
Z_true = zonotope([Z_cent,Z_gens_true]);

% check results
res = all(all(Z_del.Z == Z_true.Z));
% false if generators are not ordered in the same way ...

if res
    disp('test_deleteAligned successful');
else
    disp('test_deleteAligned failed');
end

%------------- END OF CODE --------------