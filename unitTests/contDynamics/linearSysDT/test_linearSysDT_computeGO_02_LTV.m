function res = test_linearSysDT_computeGO_02_LTV
% test_linearSysDT_computeGO_02_LTV - unit test for copmaring the GO model of an
%   LTI and an LTV system
%
% Syntax:
%    res = test_linearSysDT_computeGO_02_LTV
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       31-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10; % number of time steps

sys_LTI = linearSysDT('exampleSys',rand(3,3), rand(3,2), [], rand(2,3), rand(2,2),0.1);
p_GO_LTI = computeGO(sys_LTI, [], [], n_k);

A = cell(n_k,1);
[A{:}] = deal(sys_LTI.A);
B = cell(n_k,1);
[B{:}] = deal(sys_LTI.B);
C = cell(n_k,1);
[C{:}] = deal(sys_LTI.C);
D = cell(n_k,1);
[D{:}] = deal(sys_LTI.D);
sys_LTV = linearSysDT('exampleSys',A, B, [], C, D,0.1);
p_GO_LTV = computeGO(sys_LTV, [], [], n_k);

% check equality
res = true;
for k=1:n_k
    assertLoop(sum(abs(p_GO_LTI.A{k} - p_GO_LTV.A{k}), 'all') <= 1e-6,k)
    assertLoop(sum(abs(p_GO_LTI.C{k} - p_GO_LTV.C{k}), 'all') <= 1e-6,k)

    for j=1:k
        assertLoop(sum(abs(p_GO_LTI.B{k,j} - p_GO_LTV.B{k,j}), 'all') <= 1e-6,k,j)
        assertLoop(sum(abs(p_GO_LTI.D{k,j} - p_GO_LTV.D{k,j}), 'all') <= 1e-6,k,j)
    end
end

end

% ------------------------------ END OF CODE ------------------------------
