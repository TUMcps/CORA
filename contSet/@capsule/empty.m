function C = empty(n)
% empty - instantiates an empty capsule
%
% Syntax:
%    C = capsule.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    C - empty capsule
%
% Example: 
%    C = capsule.empty(2);
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

C = capsule(zeros(n,0),zeros(n,0),zeros(0,0));

% ------------------------------ END OF CODE ------------------------------
