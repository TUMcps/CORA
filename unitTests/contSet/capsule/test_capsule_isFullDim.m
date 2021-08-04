function res = test_capsule_isFullDim
% test_capsule_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = test_capsule_isFullDim
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% define properties
c = [2; 0; -1];
g = [-1; 1; 2];
g0 = [0; 0; 0];
r = 0.5;
r0 = 0;

% generator and radius all-zero
C = capsule(c,g0,r0);
res(1) = ~isFullDim(C);

% generator all-zero
C = capsule(c,g0,r);
res(2) = isFullDim(C);

% radius is zero
C = capsule(c,g,r0);
res(3) = ~isFullDim(C);

% generator and radius non-zero
C = capsule(c,g,r);
res(4) = isFullDim(C);


% combine results
res = all(res);

if res
    disp('test_isFullDim successful');
else
    disp('test_isFullDim failed');
end

%------------- END OF CODE --------------