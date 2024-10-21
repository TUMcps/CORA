function E_out = copy(E)
% copy - copies the ellipsoid object (used for dynamic dispatch)
%
% Syntax:
%    E_out = copy(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    E_out - copied ellipsoid object
%
% Example: 
%    E = ellipsoid(eye(2),[1;-1]);
%    E_out = copy(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
E_out = ellipsoid(E);

% ------------------------------ END OF CODE ------------------------------
