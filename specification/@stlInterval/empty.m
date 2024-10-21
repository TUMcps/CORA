function I = empty(n)
% empty - instantiates an empty STL interval
%
% Syntax:
%    I = stlInterval.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    I - empty STL interval
%
% Example: 
%    I = stlInterval.empty(1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    n = 1;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});
if n ~= 1
    throw(CORAerror('CORA:wrongValue','first','STL intervals are always one-dimensional.'));
end

I = stlInterval();

% ------------------------------ END OF CODE ------------------------------
