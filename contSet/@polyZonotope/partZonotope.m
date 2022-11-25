function [PZ,pZ_gi_cell] = partZonotope(pZ,id)
% partZonotope - computes a zonotope overapproximation in id
%
% Syntax:  
%    [PZ,pZ_gi_cell] = partZonotope(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    PZ         - zonotope overapproximation of pZ in id (in general a polyZonotope
%                 object)
%    pZ_gi_cell - cell array of polyZonotopes (generators of new zonotope
%                 factors)
%
% Example: 
%    pZ = polyZonotope(1,[1,1],[],[2,0;2,1;0,1],[1;2;3]);
%    partZono = partZonotope(pZ,2);
%    id_new = partZono.id(~ismember(partZono.id,pZ.id));
%    print(partZono,{[1;3],id_new},{'b','w'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, mptPolytope

% Author:       Victor Gassmann
% Written:      24-March-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% overapproximates pZ as a zonotope only in id_part
% pZ_gi_cell: generators of {pZ}_w which are non-constant in p
if length(id)>length(pZ.id) || ~all(ismember(id,pZ.id))
    error('provided ids not valid');
end
idv = id;
idp = pZ.id(~ismember(pZ.id,idv));
G = pZ.G;
[Id,ind] = sort(pZ.id);
expMat = pZ.expMat(ind,:);
Idv = sort(idv);
[Idp,iip] = sort(idp);
iip_rev(iip) = 1:length(iip);
ind_v = ismember(Id,Idv);
ind_p = ismember(Id,Idp);
nx = length(pZ.c);
np = length(idp);
%% compute overapproximation
% compute center
ind_nov = sum(expMat(ind_v,:),1)==0;
if all(~ind_nov)
    pZ_cp = polyZonotope(pZ.c,zeros(nx,0),zeros(nx,0),zeros(np,0),idp);
else
    eM_cp = expMat(ind_p,ind_nov); 
    pZ_cp = polyZonotope(pZ.c,G(:,ind_nov),zeros(nx,0),eM_cp(iip_rev,:),idp);
    pZ_cp = removeRedundancies(pZ_cp);
end
% compute generators
eM_vp = expMat(:,~ind_nov);
G_vp = G(:,~ind_nov);
eM_v = eM_vp(ind_v,:);
[~,~,Indc] = unique(eM_v','rows','stable');
indexE = accumarray(Indc,(1:sum(~ind_nov))',[],@(x){x});
% filter out all v-exponent vectors where the collected exponent matrix is
% all-0
indZero = cellfun(@(cc)all(all(eM_vp(ind_p,cc)==0)),indexE);
indexZero = indexE(indZero);
ii_zero = horzcat(indexZero{:});
Gzero = G_vp(:,ii_zero);
ind_even = all(mod(eM_v(:,ii_zero),2)==0,1);
c = sum(1/2*Gzero(:,ind_even),2);
Gzero(:,ind_even) = 1/2*Gzero(:,ind_even);
Ido = (max(Idp)+1:max(Idp)+length(indexE))';
I = eye(length(Ido));
PZ = polyZonotope(c,Gzero,zeros(nx,0),I(indZero,indZero),Ido(indZero));
% collect 
iE = indexE(~indZero);
pZ_gi_cell = cell(length(iE),1);
for i=1:length(iE)
    % all exponents with same v part
    eMi = eM_vp(ind_p,iE{i});
    % check if eM_vp(:,indE{i}) only dependent on v
    if all(all(eMi==0))
        tmp = polyZonotope(zeros(nx,1),sum(G_vp(:,iE{i}),2),zeros(nx,0),I(:,i),Ido);
        pZ_gi_cell{i} = polyZonotope(sum(G_vp(:,iE{i}),2),zeros(nx,0),zeros(nx,0),zeros(np,0),idp);
        PZ = exactPlus(PZ,tmp);
        continue;
    end
    % unique v exponent vector
    e_vi = eM_vp(ind_v,iE{i}(1));
    Gi = G_vp(:,iE{i});
    [eMi,Gi] = removeRedundancies(eMi,Gi);
    [eMi,Gi,ci] = removeZeroExponents(eMi,Gi);
    if sum(mod(e_vi,2))==0
        % only even exponents
        pZtmp = polyZonotope(1/2*ci,1/2*Gi,zeros(nx,0),eMi(iip_rev,:),idp);
        pZ_cp = exactPlus(pZ_cp,pZtmp);
        ci = 1/2*ci;
        Gi = 1/2*Gi;
    end
    tmp = polyZonotope(zeros(nx,1),[ci,Gi],zeros(nx,0),[zeros(np,1),eMi(iip_rev,:);repmat(I(:,i),1,1+size(eMi,2))],[idp;Ido]);
    PZ = exactPlus(PZ,tmp);
    % check if only constant ci remains
    if ~isempty(eMi)%||isempty(Gi)
        pZ_gi = polyZonotope(ci,Gi,zeros(nx,0),eMi(iip_rev,:),idp);
        pZ_gi_cell{i} = pZ_gi; 
    else
        pZ_gi_cell{i} = polyZonotope(ci,zeros(nx,0),zeros(nx,0),zeros(np,0),idp);
    end
end
% add center
PZ = exactPlus(PZ,pZ_cp);
pZ_cp_cell = {pZ_cp};
% check if pZ_cp==0
if isZero(pZ_cp)
    pZ_cp_cell = [];
end
% order important (center first important!)
pZ_gi_cell = [pZ_cp_cell;pZ_gi_cell];
%------------- END OF CODE --------------