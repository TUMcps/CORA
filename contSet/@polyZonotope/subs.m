function PZ = subs(pZ,pZin,idu)
% subs - computes substitution of pZin into pZ for elements of id idu
%
% Syntax:  
%    PZ = subs(pZ,pZin)
%    PZ = subs(pZ,pZin,idu)
%
% Inputs:
%    pZ,pZin - polyZonotope objects
%    idu - ids in pZ to be substituted
%
% Outputs:
%    PZ - resulting set as a polyZonotope object
%
% Example: 
%    % substitution
%    pZ = polyZonotope([0;1],eye(2),zeros(2,0),[2,0;0,1;1,0],[1;2;3]);
%    pZin = polyZonotope([0;0],[1,1,0;0,0,1],zeros(2,0),[1,0,0;1,0,2;0,1,0],[-1;1;2]);
%    res = subs(pZ,pZin,[1;2]);
%    print(res,{-1,[1;2],3},{'a','b','c'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: resolve

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('idu','var')
    idu = pZ.id;
end
nu = size(pZin.c,1);
nbin = size(pZin.expMat,1);
if length(idu)~=nu || ~all(ismember(idu,pZ.id))
    error('Number of outputs of pZin must match number of dep. factors specified by id.');
end
if ~isempty(pZin.Grest) && ~all(pZin.Grest==0,'all')
    error('Rest matrix of pZin is not empty or zero');
end
if isempty(pZ.G) || isempty(pZ.expMat)
    PZ = pZ;
    return;
end
%% sorting
[Id,ind] = sort(pZ.id);
[Idu,iuin] = sort(idu);
Gin = pZin.G(iuin,:);
cin = pZin.c(iuin,:);
[Idin,indin] = sort(pZin.id);
eM = pZ.expMat(ind,:);
eMin = pZin.expMat(indin,:);
% remaining ids in pZ after substitution
Idr = setdiff(Id,Idu);
% common ids in pZ.id and pZin.id that remain after substitution
Idx = Idr(ismember(Idr,Idin));
% ids not substituted but not in pZin
Idz = setdiff(Idr,Idx);
% ids in pZin but not in pZ
Idp = setdiff(Idin,Idx);
%idnew = [Idx;Idz;Idp];
%% compute substitution result
iinx = ismember(Idin,Idx);
iinp = ~iinx;
ix = ismember(Id,Idx);
iz = ismember(Id,Idz);
nu = length(Idu);
n_mons = sum(Gin~=0,2);
% check if pZin only 1 monomial per dimension
if all(cin==0) && all(n_mons<=1)
    m = size(eM,2);
    i_u = ismember(Id,Idu);
    eM_exp = eM(i_u,:);
    ind_mons = n_mons==1;
    i_emin = (Gin~=0)*(1:size(Gin,2))';
    G_subs = sum(Gin,2);
    G = pZ.G.*prod(G_subs.^eM_exp,1);
    % exponentiate
    eM_subs = zeros(length(Idin),nu);
    eM_subs(:,ind_mons) = eMin(:,i_emin(ind_mons));
    eMtmp = sum(reshape(reshape(eM_exp',[1,nu*m]).*repelem(eM_subs,1,m),...
                                [length(Idin),m,nu]),3);
    ExpMat = [eMtmp(iinx,:) + eM(ix,:);
              eM(iz,:);
              eMtmp(iinp,:)];
    
else
    degs = max(eM(ismember(Id,Idu),:),[],2);
    d = max(degs);
    Y = cell(nu,d+1);
    Gin0 = [cin,Gin];
    eMin0 = [zeros(nbin,1),eMin];
    % ^0 = 1
    Y(:,1) = repmat({{1,zeros(nbin,1)}},nu,1);
    % ^1 = original
    Y(:,2) = arrayfun(@(ii) {{Gin0(ii,:),eMin0}},(1:nu)');

    % compute powers
    for i=1:nu
        for j=3:degs(i)+1
            Y{i,j} = cell(1,2);
            [Y{i,j}{:}] = multiply(Y{i,2}{:},Y{i,j-1}{:});
        end
    end
    % multiply old monomials with new eq from pZin substituted in
    G = [];
    ExpMat = [];

    for i=1:size(pZ.G,2)
        Xprod = Y{1,eM(ismember(Id,Idu(1)),i)+1}; 
        for j=2:length(Idu)
            [Xprod{:}] = multiply(Xprod{:},Y{j,eM(ismember(Id,Idu(j)),i)+1}{:});
        end
        G = [G,pZ.G(:,i).*Xprod{1}];
        eMtmp = Xprod{2};
        eMtmp(iinx,:) = eMtmp(iinx,:) + eM(ix,i);
        eMtmp_new = [eMtmp(iinx,:);
                     repmat(eM(iz,i),1,size(eMtmp,2));
                     eMtmp(iinp,:)];
        ExpMat = [ExpMat,eMtmp_new];
    end
end
[ExpMat,G] = removeRedundantExponents(ExpMat,G);
[ExpMat,G,c] = removeZeroExponents(ExpMat,G);
%% reverse order to original
% make sure that id order is reversed again: otherwise if s.o. has some ids
% stored that are created before this function is executed, they might not 
% refer to the same dep. factors anymore
irx = ismember(Idr,Idx);
irz = ismember(Idr,Idz);
% expMat_r.id=[Idr];
expMat_r = zeros(length(Idr),size(ExpMat,2));
expMat_r(irx,:) = ExpMat(1:length(Idx),:);
expMat_r(irz,:) = ExpMat(length(Idx)+1:length(Idx)+length(Idz),:);
% expMat_p.id=[Idp];
expMat_p = ExpMat(length(Idx)+length(Idz)+1:end,:);
% now revert sort order
idr = pZ.id(~ismember(pZ.id,idu));
idp = pZin.id(ismember(pZin.id,Idp));
[~,tmp_r] = sort(idr);
[~,tmp_p] = sort(idp);
% want to go from sorted to original
indr_rev(tmp_r) = 1:length(tmp_r);
indp_rev(tmp_p) = 1:length(tmp_p);
expMat = [expMat_r(indr_rev,:);expMat_p(indp_rev,:)];
PZ = polyZonotope(c+pZ.c,G,pZ.Grest,expMat,[idr;idp]);
%------------- END OF CODE --------------