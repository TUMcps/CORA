function res = test_interval_diag
% test_interval_diag - unit test function of interval/diag
%
% Syntax:  
%    res = test_interval_diag
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
% See also: interval/diag

% Author:       Tobias Ladner
% Written:      18-July-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = [];

% 1. Empty case
I = interval();
resvec(end+1) = isequal(diag(I),I);
resvec(end+1) = isequal(diag(I,0),I);

    
% init random interval
lb = [-3; -2; -5];
ub = [4; 2; 1];
I = interval(lb,ub);

% test diagonal
D = diag(I);
resvec(end+1) = isequal(D,interval(diag(lb),diag(ub)));
resvec(end+1) = dim(I) == 3;

% test getter
resvec(end+1) = isequal(diag(D,0),I);
resvec(end+1) = isequal(diag(D,1),interval([0;0]));

% combine results
res = all(resvec);

%------------- END OF CODE --------------
