function T = empty(n)
% empty - instantiates an empty taylor model
%
% Syntax:
%    T = taylm.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    T - empty taylor model
%
% Example: 
%    T = taylm.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       13-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    n = 0;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});

if n == 0
    T = taylm(interval.empty());
    return
end

T_part = taylm(interval.empty(1));
T = T_part;
for i=2:n
    T = [T;T_part];
end

% ------------------------------ END OF CODE ------------------------------
