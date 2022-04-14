function res = test_zonotope_sampleBox
% test_zonotope_sampleBox - unit test function of sampleBox
%
% Syntax:  
%    res = test_zonotope_sampleBox
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

% Author:       Victor Gassmann
% Written:      16-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
for i=2:5:30
    G = randn(i,i);
    if rank(G)<i
        continue;
    end
    Z = zonotope([zeros(i,1),G]);
    X = sampleBox(Z,100);
    %check if all samples are contained in Z
    if ~in(Z,X)
        res = false;
        break;
    end
end
if res
    disp('test_zonotope_sampleBox successful');
else
    disp('test_zonotope_sampleBox failed');
end

%------------- END OF CODE --------------
