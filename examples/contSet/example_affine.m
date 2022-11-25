function completed = example_affine()
% example_affine - example instantiation of affine objects
%
% Syntax:  
%    completed = example_affine()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        Matthias Althoff
% Written:       29-June-2018
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% create affine object
int = interval(-1,1);
aff = affine(int);

% create taylor model object (for comparison)
maxOrder = 1;
tay = taylm(int,maxOrder,'x');

% define function
f = @(x) sin(x) * (x+1);

% evaluate the function with affine arithmetic and Taylor model
intAff = interval(f(aff))
intTay = interval(f(tay))

completed = 1;

%------------- END OF CODE --------------