function completed = example_zoo()
% example_zoo - example instantiation of affine objects
%
% Syntax:
%    completed = example_zoo()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       29-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create zoo object
int = interval(-1,1);
methods = {'interval','taylm(int)'};
maxOrder = 3;
z = zoo(int,methods,maxOrder);

% create taylor model object (for comparison)
maxOrder = 10;
tay = taylm(int,maxOrder,'x');

% define function
f = @(x) sin(x) * (x+1);

% evaluate the function with zoo-object and Taylor model
intZoo = interval(f(z))
intTay = interval(f(tay))

completed = true;

% ------------------------------ END OF CODE ------------------------------
