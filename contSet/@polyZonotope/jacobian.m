function pZdiff_cell = jacobian(pZ,id_diff)
% jacobian - computes the derivatives of each dimension of pZ with respect
% to id_diff
%
% Syntax:  
%    PZ = jacobian(pZ)
%    PZ = jacobian(pZ,id_diff)
%
% Inputs:
%    pZ      - polyZonotope object (n-dimensional)
%    id_diff - ids with respect to which pZ is differentiated
%    (d-dimensional)
%
% Outputs:
%    PZ - d-dimensional cell array where each element is a polyZonotope of
%    dimension n
%
% Example: 
%    pZ = polyZonotope(0,[1,2;3,4],[],[1,1;3,2;0,1],[1;2;3]);
%    pZ_jac_cell = jacobian(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cartProd

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('id_diff','var')
    id_diff = pZ.id;
elseif length(pZ.id)<length(id_diff) || ~all(ismember(id_diff,pZ.id))
    error('Some ids not contained in pZ.id');
end
n = length(pZ.c);
d = length(id_diff);
I = eye(length(pZ.id));
pZdiff_cell = cell(d,1);
[~,~,ii_diff] = intersect(id_diff,pZ.id,'stable');
for i=1:d
    eMtmp = pZ.expMat - I(:,ii_diff(i));
    indtmp = any(eMtmp<0,1);
    eMtmp(:,indtmp) = [];
    if ~isempty(eMtmp)
        Gtmp = pZ.expMat(ii_diff(i),~indtmp).*pZ.G(:,~indtmp);
        [E,G,c] = removeZeroExponents(eMtmp,Gtmp);
    else
        E = zeros(length(pZ.id),0);
        G = zeros(n,0);
        c = zeros(n,1);
    end
    pZdiff_cell{i} = polyZonotope(c,G,zeros(n,0),E,pZ.id);
end