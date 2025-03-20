function I = interval(O,varargin)
% interval - conversion to interval objects
%
% Syntax:
%    I = interval(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    I - interval object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/empty

% Authors:       Tobias Ladner
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval.empty(dim(O));

% ------------------------------ END OF CODE ------------------------------
