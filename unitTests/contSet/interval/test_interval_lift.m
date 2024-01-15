function res = test_interval_lift
% test_interval_lift - unit test function of lift
%
% Syntax:
%    res = test_interval_lift
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% check default case
I = interval([2;3],[4;5]);
R = lift(I,5,[1,4]);
Rtrue = interval([2;-inf;-inf;3;-inf],[4;inf;inf;5;inf]);
resvec(end+1) = isequal(R,Rtrue);

% gather results
res = all(resvec);

% empty not supported
% check empty
I = interval.empty(1);
try
    R = lift(I,2,[]);
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
