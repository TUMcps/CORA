function res = test_signal_at
% test_signal_at - unit test of signal at function
%
% Syntax:
%    res = test_signal_at
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       09-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
sig1 = signal([1.2 2.3 3.0], [true false true]);
sig2 = signal(2.7, false);

% test
res = [];
res(end+1,1) = isequal(true, sig1.at(0));
res(end+1,1) = isequal(true, sig1.at(1));
res(end+1,1) = isequal(false, sig1.at(1.2));
res(end+1,1) = isequal(false, sig1.at(2.2));
res(end+1,1) = isequal(true, sig1.at(2.3));
res(end+1,1) = isequal(true, sig1.at(3.0));

res(end+1,1) = isequal(false, sig2.at(0));
res(end+1,1) = isequal(false, sig2.at(1.9));
res(end+1,1) = isequal(false, sig2.at(2.7));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
