function res = test_taylm_interval
% test_taylm_interval - unit test of interval function
%
% Syntax:
%    res = test_taylm_interval
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none

% Authors:       Tobias Ladner
% Written:       17-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init Taylor model
tx = taylm(interval(1,4),4,'x');
ty = taylm(interval(1,4),4,'y');

% evaluate on sin
tay = sin(tx + ty);

% test all options
I = interval(tay); % default 'int'
resvec(end+1) = isequal(I,interval(-10.5053,9.1269),1e-4);
I = interval(tay,'int');
resvec(end+1) = isequal(I,interval(-10.5053,9.1269),1e-4);
I = interval(tay,'bnb');
resvec(end+1) = isequal(I,interval(-5.5094,4.3090),1e-4);
I = interval(tay,'bnbAdv');
resvec(end+1) = isequal(I,interval(-3.3641,2.8901),1e-4);
I = interval(tay,'bernstein');
resvec(end+1) = isequal(I,interval(-4.2425,3.5428),1e-4);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
