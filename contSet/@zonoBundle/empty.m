function zB = empty(n)
% empty - instantiates an empty zonotope bundle
%
% Syntax:
%    zB = zonoBundle.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    zB - empty zonotope bundle
%
% Example: 
%    zB = zonoBundle.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

zB = zonoBundle({zonotope(zeros(n,0))});

% ------------------------------ END OF CODE ------------------------------