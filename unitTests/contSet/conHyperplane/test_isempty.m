function res = test_isempty
% test_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_isempty
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  03-May-2020 (add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

hyp1 = conHyperplane(halfspace([0 1],3),eye(2),5);
hyp2 = conHyperplane(halfspace([0 1],3),eye(2),0);
hyp3 = conHyperplane();

% combine tests
res = ~isempty(hyp1) && ~isempty(hyp2) && isempty(hyp3);

if res
    disp('test_isempty successful');
else
    disp('test_isempty failed');
end

%------------- END OF CODE --------------