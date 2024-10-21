function res = isBounded(SpS)
% isBounded - determines whether a spectrahedron is bounded and saves it in
%    the hidden propertey "ESumRep"
%
% Syntax:
%    res = isBounded(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    res - true/false whether the spectrahedral shadow is bounded
%
% Example: 
%    SpS_bounded = spectraShadow([1 0 -1 0;0 -1 0 1]);
%    res = isBounded(SpS_bounded)
%    SpS_unbounded = spectraShadow([zeros(3) zeros(3)]);
%    res = isBounded(SpS_unbounded)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       06-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Deal with trivial case
if isemptyobject(SpS)
    res = true;
    return
end

if ~isempty(SpS.bounded.val)
    res = SpS.bounded.val;
    return
end

% Assume it is bounded
res = true;

% Loop over dimensions, check how far we can go in each axis direction,
% with the support function
n = dim(SpS);
for i=1:n
    ei = zeros([n 1]);
    ei(i) = 1;
    val_p = supportFunc_(SpS,ei,'upper');
    val_m = supportFunc_(SpS,-ei,'upper');
    
    if (val_p == Inf) || (val_m == Inf)
        res = false;
        SpS.bounded.val = false;
        return
    end
end

SpS.bounded.val = true;


% ------------------------------ END OF CODE ------------------------------
