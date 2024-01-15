function ls = Inf(n)
% Inf - instantiates a fullspace level set
%
% Syntax:
%    ls = levelSet.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    ls - fullspace level set
%
% Example: 
%    ls = levelSet.Inf(2);
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

vars = sym('x',[n,1]);
eq = 0*vars - 1;
ls = levelSet(eq,vars,{"<="});

% ------------------------------ END OF CODE ------------------------------
