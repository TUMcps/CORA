function SpS_out = enclose(SpS,varargin)
% enclose - encloses a spectrahedral shadow and its affine transformation
%
% Description:
%    Computes the set
%    {a x1 + (1 - a) * x2 | x1 \in SpS, x2 \in SpS2, a \in [0,1]}
%    where SpS2 = M*SpS + SpSplus
%
% Syntax:
%    SpS_out = enclose(SpS,SpS2)
%    SpS_out = enclose(SpS,M,SpSplus)
%
% Inputs:
%    SpS - spectraShadow object
%    SpS2 - spectraShadow object
%    M - matrix for the linear transformation
%    SpSplus - spectraShadow object added to the linear transformation
%
% Outputs:
%    SpS_out - spectraShadow object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS1 = spectraShadow([A0 A1 A2]);
%    M = [-1 0; 0 -1];
%    SpS2 = M*SpS1 + [0.5;0.5];
%    SpS=enclose(SpS1,SpS2);
%
%    figure; hold on;
%    plot(SpS1,[1,2],'r');
%    plot(SpS2,[1,2],'g');
%    plot(SpS,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/enclose, zonotope/enclose

% Authors:       Adrian Kulmburg
% Written:       17-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
if nargin == 2
    SpS2 = varargin{1};
elseif nargin == 3
    SpS2 = (varargin{1}*SpS) + varargin{2};
end

SpS_out = convHull(SpS, SpS2);

% ------------------------------ END OF CODE ------------------------------
