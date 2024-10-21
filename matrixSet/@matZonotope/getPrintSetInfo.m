function [abbrev,propertyOrder] = getPrintSetInfo(S)
% getPrintSetInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSetInfo(S)
%
% Inputs:
%    S - contSet or matrixSet
%
% Outputs:
%    abbrev - set abbreviation
%    propertyOrder - order of the properties
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: printSet

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

abbrev = 'matZ';
propertyOrder = {'C','G'};

% ------------------------------ END OF CODE ------------------------------
