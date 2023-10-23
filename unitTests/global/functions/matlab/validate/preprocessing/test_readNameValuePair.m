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
res = val == 12;
res(end+1,1) = length(NVpairs_) == 4;

% read value for name 'Name', check case insensitivity
[NVpairs_,val] = readNameValuePair(NVpairs,'NAME');
res(end+1,1) = strcmp(val,'John');
res(end+1,1) = length(NVpairs_) == 4;

% read value for name 'Children', spelt wrongly
[NVpairs_,val] = readNameValuePair(NVpairs,'Childrn');
res(end+1,1) = isempty(val);
res(end+1,1) = length(NVpairs_) == 6;

% read out existing name with default value
[NVpairs_,val] = readNameValuePair(NVpairs,'Age','isscalar',2);
res(end+1,1) = val == 12;
res(end+1,1) = length(NVpairs_) == 4;

% read out non-existing name with default value
[NVpairs_,val] = readNameValuePair(NVpairs,'Parents','isscalar',2);
res(end+1,1) = val == 2;
res(end+1,1) = length(NVpairs_) == 6;

% fail the check
try 
    [NVpairs_,val] = readNameValuePair(NVpairs,'Age','ischar');
    res(end+1,1) = false;
end
try 
    [NVpairs_,val] = readNameValuePair(NVpairs,'Name',@isscalar);
    res(end+1,1) = false;
end


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
