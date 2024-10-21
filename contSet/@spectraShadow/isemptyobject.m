function res = isemptyobject(SpS)
% isemptyobject - checks if a spectrahedral shadow is fully empty
%
% Syntax:
%    res = isemptyobject(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    res - true/false
%
% Example: 
%    SpS = spectraShadow([eye(3) eye(3)]);
%    isemptyobject(SpS) % false
%    SpS = spectraShadow.empty()
%    isemptyobject(SpS) % true
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no center vector and generators
% (by construction, SpS.A can never be empty)
res = isempty(SpS.c) && isempty(SpS.G);

% ------------------------------ END OF CODE ------------------------------
