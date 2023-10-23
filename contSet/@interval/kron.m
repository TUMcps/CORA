function res = kron(I1,I2)
% kron - Kronecker product of two intervals
%
% Syntax:
%    intMatKron = kron(I1,I2)
%
% Inputs:
%    I1,I2 - intervals
%
% Outputs:
%    intMatKron - interval matrix 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: times

% Authors:       Victor Gassmann
% Written:       22-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Vals = [];
if ~isa(I1,'interval')
    sz = size(I1).*size(I2.inf);
    pv1_mat = kron(I1,I2.inf); Vals(:,1) = pv1_mat(:);
    pv2_mat = kron(I1,I2.sup); Vals(:,2) = pv2_mat(:);
elseif ~isa(I2,'interval')
    sz = size(I1.inf).*size(I2);
    pv1_mat = kron(I1.inf,I2); Vals(:,1) = pv1_mat(:);
    pv2_mat = kron(I1.sup,I2); Vals(:,2) = pv2_mat(:);
else
    sz = size(I1.inf).*size(I2.inf);
    % possible combinations
    pv1_mat = kron(I1.inf,I2.inf); Vals(:,1) = pv1_mat(:);
    pv2_mat = kron(I1.inf,I2.sup); Vals(:,2) = pv2_mat(:);
    pv3_mat = kron(I1.sup,I2.inf); Vals(:,3) = pv3_mat(:);
    pv4_mat = kron(I1.sup,I2.sup); Vals(:,4) = pv4_mat(:);
end

% to find min and max
inf_mat = reshape(min(Vals,[],2),sz);
sup_mat = reshape(max(Vals,[],2),sz);

res = interval(inf_mat,sup_mat);

% ------------------------------ END OF CODE ------------------------------
