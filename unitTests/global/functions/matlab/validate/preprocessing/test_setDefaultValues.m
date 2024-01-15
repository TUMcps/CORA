function res = test_setDefaultValues
% test_setDefaultValues - unit test function for setting of default values
%
% Syntax:
%    res = test_setDefaultValues()
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

% some user-provided values and default values
givenValues = {10,'bbb',polyZonotope.empty(2)};
defaultValues = {5,'aaa',zonotope.empty(2)};

% no input arguments given
[val1,val2,val3] = setDefaultValues(defaultValues,{});
res = val1 == defaultValues{1};
res(end+1,1) = strcmp(val2,defaultValues{2});
res(end+1,1) = isequal(val3,defaultValues{3});

% one input argument given
[val1,val2,val3] = setDefaultValues(defaultValues,givenValues(1));
res(end+1,1) = val1 == givenValues{1};
res(end+1,1) = strcmp(val2,defaultValues{2});
res(end+1,1) = isequal(val3,defaultValues{3});

% all input arguments given
[val1,val2,val3] = setDefaultValues(defaultValues,givenValues);
res(end+1,1) = val1 == givenValues{1};
res(end+1,1) = strcmp(val2,givenValues{2});
res(end+1,1) = isequal(val3,givenValues{3});

% too many input arguments given
[val1,val2,val3] = setDefaultValues(defaultValues,[givenValues,givenValues]);
res(end+1,1) = val1 == givenValues{1};
res(end+1,1) = strcmp(val2,givenValues{2});
res(end+1,1) = isequal(val3,givenValues{3});

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
