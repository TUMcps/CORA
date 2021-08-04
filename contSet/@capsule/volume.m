function vol = volume(C)
% volume - Computes the volume of a capsule
%
% Syntax:  
%    vol = volume(C)
%
% Inputs:
%    C - capsule
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
% See also: none

% Author:       Matthias Althoff
% Written:      05-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check empty set
if isempty(C)
    vol = []; return;
end

% the volume is obtained by the volume of an n-dimensional ball and an
% n-dimensional cylinder

% dimension
d = length(C.c);

% volume of n-dimenional ball; requires Leonhard Euler's gamma function
volBall = pi^(d/2)/gamma(d/2+1)*C.r^d;

% volume of n-dimensional cylinder is volume of (n-1)-dimensional ball
% times length
volBall_minus1 = pi^((d-1)/2)/gamma((d-1)/2+1)*C.r^(d-1);
volCylinder = volBall_minus1*(2*norm(C.g));

% total volume
vol = volBall + volCylinder;


%------------- END OF CODE --------------