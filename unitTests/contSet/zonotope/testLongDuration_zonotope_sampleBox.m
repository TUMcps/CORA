function res = testLongDuration_zonotope_sampleBox
% testLongDuration_zonotope_sampleBox - unit test function of sampleBox
%
% Syntax:  
%    res = testLongDuration_zonotope_sampleBox
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

% Author:       Victor Gassmann
% Written:      16-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% loop over varying dimension
for i=2:5:30
    
    % random generator matrix
    G = randn(i,i);

    % check rank (should not be smaller than i)
    if rank(G) < i
        continue;
    end

    % instantiate zonotope
    Z = zonotope([zeros(i,1),G]);

    % sample
    X = sampleBox(Z,100);

    % check if all samples are contained in Z
    if ~contains(Z,X)
        res = false; return
    end
end

%------------- END OF CODE --------------