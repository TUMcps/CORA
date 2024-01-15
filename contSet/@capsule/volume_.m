function vol = volume_(C,varargin)
% volume_ - Computes the volume of a capsule
%
% Syntax:
%    vol = volume_(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    vol - volume
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    vol = volume(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Matthias Althoff
% Written:       05-March-2019
% Last update:   18-August-2022 (MW, include standardized preprocessing)
% Last revision: 27-March-2023 (MW, rename volume_)

% ------------------------------ BEGIN CODE -------------------------------

% the volume is obtained by the volume of an n-dimensional ball and an
% n-dimensional cylinder

% dimension
n = dim(C);

if representsa_(C,'emptySet',eps)
    vol = 0; return
end

% special case
if n == 0
    vol = 0; return
end

% volume of n-dimenional ball; requires Euler's gamma function
volBall = pi^(n/2)/gamma(n/2+1)*C.r^n;

% volume of n-dimensional cylinder is volume of (n-1)-dimensional ball
% times length
volBall_minus1 = pi^((n-1)/2)/gamma((n-1)/2+1)*C.r^(n-1);
volCylinder = volBall_minus1*(2*norm(C.g));

% total volume
vol = volBall + volCylinder;

% ------------------------------ END OF CODE ------------------------------
