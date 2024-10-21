function [abbrev,propertyOrder] = getPrintSetInfo(P)
% getPrintSetInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSetInfo(P)
%
% Inputs:
%    P - polytope
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

abbrev = 'P';
if P.isHRep.val
    propertyOrder = {'A','b','Ae','be'};
else % isVrep
    propertyOrder = {'V'};
end

% ------------------------------ END OF CODE ------------------------------
