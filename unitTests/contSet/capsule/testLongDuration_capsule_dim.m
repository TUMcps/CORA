function res = testLongDuration_capsule_dim
% testLongDuration_capsule_dim - unit test function of dim
%
% Syntax:  
%    res = testLongDuration_capsule_dim
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
% Written:      27-Sep-2019
% Last update:  12-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

% Random
C = cell(10,1);
dimension = zeros(10,1);
for i=1:10
    dimension(i) = ceil(10*rand(1));
    center = rand(dimension(i),1);
    generator = rand(dimension(i),1);
    r = rand(1);
    C{i} = capsule(center,generator,r);
end

res = true;
for i=1:length(C)
    res = res && dim(C{i}) == dimension(i);
end

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------