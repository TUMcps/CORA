function res = test_linErrorBound_linErrorBound
% test_linErrorBound_linErrorBound - unit test for helper class
%    linErrorBound
%
% Syntax:
%    res = test_linErrorBound_linErrorBound
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       10-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% correct initialization
emax = 1; tFinal = 1;
errs = linErrorBound(emax,tFinal);
assert(errs.emax == emax);
assert(errs.tFinal == tFinal);


% wrong initializations
% - wrong number of input arguments
assertThrowsAs(@linErrorBound,'MATLAB:narginchk:notEnoughInputs');
assertThrowsAs(@linErrorBound,'MATLAB:narginchk:notEnoughInputs',1);
assertThrowsAs(@linErrorBound,'MATLAB:TooManyInputs',1,1,1);
% - wrong values
assertThrowsAs(@linErrorBound,'CORA:wrongInputInConstructor',-1,1);
assertThrowsAs(@linErrorBound,'CORA:wrongInputInConstructor',1,-1);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
