function [abbrev,propertyOrder] = getPrintSetInfo(pZ)
% getPrintSetInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSetInfo(pZ)
%
% Inputs:
%    pZ - polyZonotope
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

abbrev = 'pZ';
propertyOrder = {'c','G','GI','E','id'};

% ------------------------------ END OF CODE ------------------------------
