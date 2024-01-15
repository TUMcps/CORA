function res = test_conHyperplane_representsa
% test_conHyperplane_representsa - unit test function of representsa
%
% Syntax:
%    res = test_conHyperplane_representsa
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

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       17-September-2019
% Last update:   03-May-2020 (add empty case)
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. comparison to empty set
hyp = conHyperplane([0 1],3,eye(2),[1;1]);
res(end+1,1) = representsa(hyp,'emptySet');
hyp = conHyperplane([0 1],3,eye(2),[1;4]);
res(end+1,1) = ~representsa(hyp,'emptySet');


% 2. comparison to interval
hyp = conHyperplane([-1 0],1);
res(end+1,1) = representsa(hyp,'interval');
hyp = conHyperplane([-1 0],1,[0 1; 0 -1],[1; 1]);
res(end+1,1) = representsa(hyp,'interval');


% 3. comparison to hyperplane
% empty constraint set
a = [3; 2; -1]; b = 0.5;
C = [0 0 0]; d = 1;
hyp = conHyperplane(a,b,C,d);
res(end+1,1) = representsa(hyp,'hyperplane');

% constraint set not empty but still unconstrained hyperplane
I = eye(2); a = 0.7;
hyp = conHyperplane(I(1,:),0,[I(1,:);-I(1,:)],[a;a]);
res(end+1,1) = representsa(hyp,'hyperplane');
    
% constraint set which actually constrains hyperplane
% could be 1D and then not consistent (see "conHyperplane"
% constructor); otherwise trivially always an unconstrained hyperplane
d = [1; -1]; d = d/norm(d);
hyp = conHyperplane(d,1,[I(1,:);-I(1,:)],ones(2,1));
res(end+1,1) = ~representsa(hyp,'hyperplane');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
