function [propertyOrder] = getPrintSystemInfo(nlnARX)
% getPrintSystemInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSystemInfo(nlnARX)
%
% Inputs:
%    nlnARX - nonlinearARX object
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

% 'number of past time steps' is stored in 'nrOfStates' ...
propertyOrder = {'name','mFile','dt','nrOfOutputs','nrOfInputs','nrOfStates'};

% ------------------------------ END OF CODE ------------------------------
