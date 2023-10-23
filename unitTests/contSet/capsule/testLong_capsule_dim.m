function res = testLong_capsule_dim
% testLong_capsule_dim - unit test function of dim
%
% Syntax:
%    res = testLong_capsule_dim
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

% Authors:       Mark Wetzlinger
% Written:       27-September-2019
% Last update:   12-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random
C = cell(10,1);
n = zeros(10,1);
for i=1:10
    n(i) = ceil(10*rand(1));
    center = rand(n(i),1);
    generator = rand(n(i),1);
    r = rand(1);
    C{i} = capsule(center,generator,r);
end

res = true;
for i=1:length(C)
    res = res && dim(C{i}) == n(i);
end

% ------------------------------ END OF CODE ------------------------------
