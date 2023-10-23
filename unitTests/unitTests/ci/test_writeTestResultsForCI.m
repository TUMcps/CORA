function res = test_writeTestResultsForCI
% test_writeTestResultsForCI - unit test function for writeTestResultsForCI
% 
% 
% Syntax:
%    res = test_writeTestResultsForCI
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
% Written:       22-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

writeTestResultsForCI('short');

% tes result text
text = fileread('resultText.txt');
% not empty
resvec(end+1) = ~isempty(text);
% no new line
resvec(end+1) = length(split(text, "\n")) == 1;

% test exit code
text = fileread('failed.txt');
% just single number
resvec(end+1) = ~isnan(str2double(text));

% delete files
delete resultText.txt;
delete failed.txt;

% return result
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
