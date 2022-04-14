function res = test_polyZonotope_isequal
% test_polyZonotope_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_polyZonotope_isequal
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
% Written:      01-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create polyZonotopes
c = [2; -3];
G1 = [2, 3;
     -1, 0];
G1rest = [1, 2, 4;
          5, 6, 0];
expMat1 = [2, 0;
           0, 1];
pZ1 = polyZonotope(c,G1,G1rest,expMat1);
G2 = G1 + ones(2,2);
pZ2 = polyZonotope(c,G2,G1rest,expMat1);
G3rest = [1, 2, 0, 4;
          5, 6, 0, 0];
pZ3 = polyZonotope(c,G1,G3rest,expMat1);

% check result
res = isequal(pZ1,pZ3) && ~isequal(pZ1,pZ2);

if res
    disp('test_polyZonotope_isequal successful');
else
    disp('test_polyZonotope_isequal failed');
end

%------------- END OF CODE --------------