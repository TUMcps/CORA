function res = test_conHyperplane_isempty
% test_conHyperplane_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_conHyperplane_isempty
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

% Author:       Mark Wetzlinger
% Written:      17-September-2019
% Last update:  03-May-2020 (add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

hyp1 = conHyperplane(halfspace([0 1],3),eye(2),ones(2,1));
hyp2 = conHyperplane(halfspace([0 1],3),eye(2),zeros(2,1));
hyp3 = conHyperplane();

% combine tests
res = ~isempty(hyp1) && ~isempty(hyp2) && isempty(hyp3);

%------------- END OF CODE --------------