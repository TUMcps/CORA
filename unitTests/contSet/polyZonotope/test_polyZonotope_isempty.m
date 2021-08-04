function res = test_polyZonotope_isempty
% test_polyZonotope_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_polyZonotope_isempty
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
G = [2, 3;
     -1, 0];
Grest = [1, 2, 4;
          5, 6, 0];
expMat = [2, 0;
           0, 1];
pZ1 = polyZonotope(c,G,Grest,expMat);
pZ2 = polyZonotope(c,G,[],expMat);
pZ3 = polyZonotope();

% check result
res = ~isempty(pZ1) && ~isempty(pZ2) && isempty(pZ3);

if res
    disp('test_polyZonotope_isempty successful');
else
    disp('test_polyZonotope_isempty failed');
end

%------------- END OF CODE --------------