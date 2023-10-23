function E = ellipsoid(I)
% ellipsoid - Converts an interval object into an ellipsoid object
%
% Syntax:
%    E = ellipsoid(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    I = interval([1;-1], [2; 1]);
%    E = ellipsoid(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Victor Gassmann
% Written:       15-October-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert interval to zonotope (exact computation very efficient because
% generator matrix is square)
E = ellipsoid(zonotope(I),'outer:exact');

% ------------------------------ END OF CODE ------------------------------
