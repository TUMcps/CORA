function res = checkNameValuePairs(NVpairs,fullList)
% checkNameValuePairs - checks if every name from a given list of
%    name-value pairs is part of the full list of admissible name-value
%    pairs (check is case-insensitive and works for chars and strings)
%
% Syntax:
%    checkNameValuePairs(NVpairs,fullList)
%
% Inputs:
%    NVpairs - cell-array: list of given name-value pairs
%    fullList - cell-array: list of all potential name-value pairs
%
% Outputs:
%    res = true/false
%
% Example: 
%    NVpairs = {'Splits',8,'Order',5};
%    checkNameValuePairs(NVpairs,{'Splits','Order'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: readNameValuePairs

% Authors:       Mark Wetzlinger
% Written:       08-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default result
res = true;

% number of given name-value pairs
nrNVpairs = length(NVpairs) / 2;

% read list of name of given name-value pairs
names = NVpairs(1:2:length(NVpairs));

% go through all given names and check whether for containment in full list
for i=1:nrNVpairs
    if ~ismember(lower(names{i}),lower(fullList))
        throw(CORAerror('CORA:redundantNameValuePair',names{i},fullList));
    end
end

% ------------------------------ END OF CODE ------------------------------
