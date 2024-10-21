function [abbrev,propertyOrder] = getPrintSetInfo(SpS)
% getPrintSetInfo - returns all information to properly print a set 
%    to the command window 
%
% Syntax:
%    [abbrev,propertyOrder] = getPrintSetInfo(SpS)
%
% Inputs:
%    SpS - spectraShadow
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

abbrev = 'SpS';
propertyOrder = {'A','c','G'};

% ------------------------------ END OF CODE ------------------------------
