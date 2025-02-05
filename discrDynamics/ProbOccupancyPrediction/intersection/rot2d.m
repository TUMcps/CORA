function pts_out = rot2d(theta, pts_in)
% rot2d - rotates a set of points of frame in to frame out
%
% Syntax:
%    pts_out = rot2d(theta, pts_in)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       ???
% Last update:   19-August-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

assert(size(pts_in,1)==2);
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
pts_out = R*pts_in;

end

% ------------------------------ END OF CODE ------------------------------
