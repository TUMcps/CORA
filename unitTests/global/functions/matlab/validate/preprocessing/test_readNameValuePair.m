function res = test_readNameValuePair
% test_readNameValuePair - unit test function for reading of name-value pairs
%
% Syntax:
%    res = test_readNameValuePair()
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init name-value pairs
NVpairs = {'Age',12,'Name','John','Children',2};

% read value for name 'Age'
[NVpairs_,val] = readNameValuePair(NVpairs,'Age');
assert(val == 12);
assert(length(NVpairs_) == 4);

% read value for name 'Name', check case insensitivity
[NVpairs_,val] = readNameValuePair(NVpairs,'NAME');
assert(strcmp(val,'John'));
assert(length(NVpairs_) == 4);

% read value for name 'Children', spelt wrongly
[NVpairs_,val] = readNameValuePair(NVpairs,'Childrn');
assert(isempty(val));
assert(length(NVpairs_) == 6);

% read out existing name with default value
[NVpairs_,val] = readNameValuePair(NVpairs,'Age','isscalar',2);
assert(val == 12);
assert(length(NVpairs_) == 4);

% read out non-existing name with default value
[NVpairs_,val] = readNameValuePair(NVpairs,'Parents','isscalar',2);
assert(val == 2);
assert(length(NVpairs_) == 6);

% fail the check
assertThrowsAs(@readNameValuePair,'CORA:specialError',NVpairs,'Age','ischar');
assertThrowsAs(@readNameValuePair,'CORA:specialError',NVpairs,'Name','isscalar');


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
