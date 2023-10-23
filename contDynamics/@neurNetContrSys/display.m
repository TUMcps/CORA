function display(obj)
% display - displays the properties of a neurNetContrSys object
%
% Syntax:
%    display(obj)
%
% Inputs:
%    obj - neurNetContrSys object
%
% Outputs:
%    -
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

disp("sys:")
display(obj.sys)

disp("nn:")
display(obj.nn)
disp(" ")

disp("dt: (neural network sampling time)")
disp(obj.dt)

% ------------------------------ END OF CODE ------------------------------
