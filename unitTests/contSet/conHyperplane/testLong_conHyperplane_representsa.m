function res = testLong_conHyperplane_representsa
% testLong_conHyperplane_representsa - unit test function of representsa
%
% Syntax:
%    res = testLong_conHyperplane_representsa
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

% Authors:       Victor Gassmann
% Written:       22-March-2021
% Last update:   --- 
% Last revision: 26-July-2023 (MW, rename "...representsa")

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = true;

% number of tests
nrTests = 100;
for i=1:nrTests

    % random dimension
    n = randi(30);

    % empty constraint set
    hyp1 = conHyperplane(randn(1,n),randn,zeros(1,n),abs(randn));

    % constraint set not empty but still unconstrained hyperplane
    I = eye(n);
    a = abs(randn);
    hyp2 = conHyperplane(I(1,:),0,[I(1,:);-I(1,:)],[a;a]);
    if ~representsa(hyp1,'hyperplane') || ~representsa(hyp2,'hyperplane')
        throw(CORAerror('CORA:testFailed'));
    end

    % constraint set which actually constrains hyperplane
    % could be 1D and then not consistent (see "conHyperplane"
    % constructor); otherwise trivially always an unconstrained hyperplane
    if n > 1
        d = randn(1,n); d = d/norm(d);
        hyp3 = conHyperplane(d,1,[I;-I],ones(2*n,1));
        if representsa(hyp3,'hyperplane')
            throw(CORAerror('CORA:testFailed'));
        end
    end
    
end

% ------------------------------ END OF CODE ------------------------------
