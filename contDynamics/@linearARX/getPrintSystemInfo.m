function [propertyOrder] = getPrintSystemInfo(linARX)
% getPrintSystemInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSystemInfo(linARX)
%
% Inputs:
%    linARX - linearARX object
%
% Outputs:
%    propertyOrder - order of the properties
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/printSystem

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

propertyOrder = {'name','A_bar','B_bar','dt'};

% ------------------------------ END OF CODE ------------------------------
