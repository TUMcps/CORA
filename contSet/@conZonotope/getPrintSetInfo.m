function [abbrev,propertyOrder] = getPrintSetInfo(cZ)
% getPrintSetInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSetInfo(cZ)
%
% Inputs:
%    cZ - constrained zonotope
%
% Outputs:
%    abbrev - set abbreviation
%    propertyOrder - order of the properties
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

abbrev = 'cZ';
propertyOrder = {'c','G','A','b'};

% ------------------------------ END OF CODE ------------------------------
