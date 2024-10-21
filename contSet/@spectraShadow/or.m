function S_out = or(SpS,S,varargin)
% or - returns the convex hull over the union of a spectrahedral shadow and
%    a contSet object
%
% Syntax:
%    S_out = SpS | S
%    S_out = or(SpS,S)
%    S_out = or(SpS,S,mode)
%
% Inputs:
%    SpS - spectraShadow object
%    S - contSet object
%    mode - 'exact', 'outer', 'inner'
%
% Outputs:
%    S_out - convex hull over the union of S1 and S2
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    SpS_convHull = SpS | (SpS+[1;1]);
%    
%    figure; hold on;
%    plot(SpS);
%    plot(SpS+[1;1]);
%    plot(SpS_convHull);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call convex hull method
S_out = convHull_(SpS,S,varargin{:});
    
% ------------------------------ END OF CODE ------------------------------
