function hs = Inf(n)
% Inf - instantiates a fullspace halfspace object
%
% Syntax:
%    hs = halfspace.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    hs - fullspace halfspace object
%
% Example: 
%    hs = halfspace.Inf(2);
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

% halfspace 0*x <= 1 is fulfilled for all values for x
hs = halfspace(zeros(1,n),1);

% ------------------------------ END OF CODE ------------------------------
