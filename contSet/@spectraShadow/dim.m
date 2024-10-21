function n = dim(SpS)
% dim - returns the dimension of a spectrahedral shadow
%
% Syntax:
%    n = dim(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    n - dimension of the spectrahedron SpS
%
% Example:
%    SpS = spectraShadow([eye(2) eye(2)]);
%    dim(SpS) % Should be 1
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       24-April-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    c = SpS.c;
    n = size(c,1);
catch ME
    if isemptyobject(SpS)
        n = 0;
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
