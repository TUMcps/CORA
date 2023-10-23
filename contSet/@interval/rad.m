function r = rad(I)
% rad - returns the radius of an interval
%
% Syntax:
%    r = rad(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    r - numerical value
%
% Example: 
%    I = interval([-1 1], [1 2]);
%    r = rad(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       26-June-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

r = 0.5*(I.sup - I.inf);

% ------------------------------ END OF CODE ------------------------------
