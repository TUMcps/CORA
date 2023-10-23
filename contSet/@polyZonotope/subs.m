function pZ = subs(pZ,pZin,varargin)
% subs - computes substitution of pZin into pZ for elements of id idu
%
% Syntax:
%    pZ = subs(pZ,pZin)
%    pZ = subs(pZ,pZin,idu)
%
% Inputs:
%    pZ,pZin - polyZonotope objects
%    idu - (optional) ids in pZ to be substituted
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope([0;1],eye(2),zeros(2,0),[2,0;0,1;1,0],[1;2;3]);
%    pZin = polyZonotope([0;0],[1,1,0;0,0,1],zeros(2,0),[1,0,0;1,0,2;0,1,0],[-1;1;2]);
%    res = subs(pZ,pZin,[1;2]);
%    print(res,'Ids',{-1,[1;2],3},'Vars',{'a','b','c'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: resolve

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% input argument check
% default value
idu = pZ.id;

if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
if nargin == 3 && ~isempty(varargin{1})
    if ~isnumeric(varargin{1}) || ~isvector(varargin{1}) || any(mod(varargin{1},1) ~= 0)
        throw(CORAerror('CORA:wrongValue','third','Identifiers have to be correct'));
    end
    idu = varargin{1};
end

nu = size(pZin.c,1);
nbin = size(pZin.E,1);
if length(idu)~=nu || ~all(ismember(idu,pZ.id))
    throw(CORAerror('CORA:wrongValue','third',...
        'Number of outputs of pZin must match number of dep. factors specified by id.'));
end
if ~isempty(pZin.GI) && ~all(pZin.GI==0,'all')
    throw(CORAerror('CORA:wrongValue','second',...
        'Rest matrix of pZin is not empty or zero'));
end

%% trivial case
if isempty(pZ.G) || isempty(pZ.E)
    return;
end

%% sorting
% sort to make sure that non-sorted input idu etc does not cause problems
[Id,ind] = sort(pZ.id);
[Idu,iuin] = sort(idu);
Gin = pZin.G(iuin,:);
cin = pZin.c(iuin,:);
[Idin,indin] = sort(pZin.id);
eM = pZ.E(ind,:);
eMin = pZin.E(indin,:);
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
% get masks of ids
ind_inx = ismember(Idin,Idx);
ind_inp = ~ind_inx;
ix = ismember(Id,Idx);
iz = ismember(Id,Idz);
nu = length(Idu);
n_mons = sum(Gin~=0,2);
% check if pZin only 1 monomial per dimension
if all(cin==0) && all(n_mons<=1)
    % more efficient computation if all dimensions of pZin only contain one
    % monomial
    m = size(eM,2);
    i_u = ismember(Id,Idu);
    eM_exp = eM(i_u,:);
    ind_mons = n_mons==1;
    i_emin = (Gin~=0)*(1:size(Gin,2))';
    % always at most one entry ~=0, so sum used to get matrix to "vector"
    % form
    G_subs = sum(Gin,2);
    % compute resulting generator matrix
    G = pZ.G.*prod(G_subs.^eM_exp,1);
    % exponentiate
    eM_subs = zeros(length(Idin),nu);
    eM_subs(:,ind_mons) = eMin(:,i_emin(ind_mons));
    eMtmp = sum(reshape(reshape(eM_exp',[1,nu*m]).*repelem(eM_subs,1,m),...
                                [length(Idin),m,nu]),3);
    % combine result
    E = [eMtmp(ind_inx,:) + eM(ix,:);
              eM(iz,:);
              eMtmp(ind_inp,:)];
    
else
    % extract max degree for each idu entry
    degs = max(eM(ismember(Id,Idu),:),[],2);
    d = max(degs);
    % since we want to replace a dep. factor of pZ with the corresponding
    % 1D polyZonotope from pZin, we need all solutions of
    % project(pZin,i)^(0:d(i))
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
    E = [];
    
    for i=1:size(pZ.G,2)
        % extract the appropriate power
        Xprod = Y{1,eM(ismember(Id,Idu(1)),i)+1}; 
        for j=2:length(Idu)
            % multiply the result of each separate idu entry
            [Xprod{:}] = multiply(Xprod{:},Y{j,eM(ismember(Id,Idu(j)),i)+1}{:});
        end
        % the result of the monomial computation (so substituting the
        % dimensions of pZin into dep. factors of pZ) still
        % requires the computation with the coefficients of pZ, i.e., G
        G = [G,pZ.G(:,i).*Xprod{1}];
        % construct exponent matrix
        eMtmp = Xprod{2};
        eMtmp(ind_inx,:) = eMtmp(ind_inx,:) + eM(ix,i);
        eMtmp_new = [eMtmp(ind_inx,:);
                     repmat(eM(iz,i),1,size(eMtmp,2));
                     eMtmp(ind_inp,:)];
        E = [E,eMtmp_new];
    end
end
[E,G] = removeRedundantExponents(E,G);
[E,G,c] = removeZeroExponents(E,G);

%% reverse order to original
% make sure that id order is reversed again: otherwise if s.o. has some ids
% stored that are created before this function is executed, they might not 
% refer to the same dep. factors anymore
irx = ismember(Idr,Idx);
irz = ismember(Idr,Idz);
% E_r.id=[Idr];
E_r = zeros(length(Idr),size(E,2));
E_r(irx,:) = E(1:length(Idx),:);
E_r(irz,:) = E(length(Idx)+1:length(Idx)+length(Idz),:);
% E_p.id=[Idp];
E_p = E(length(Idx)+length(Idz)+1:end,:);
% now revert sort order
idr = pZ.id(~ismember(pZ.id,idu));
idp = pZin.id(ismember(pZin.id,Idp));
[~,tmp_r] = sort(idr);
[~,tmp_p] = sort(idp);
% want to go from sorted to original
indr_rev(tmp_r) = 1:length(tmp_r);
indp_rev(tmp_p) = 1:length(tmp_p);
E = [E_r(indr_rev,:);E_p(indp_rev,:)];
% instantiate result
pZ = polyZonotope(c+pZ.c,G,pZ.GI,E,[idr;idp]);

% ------------------------------ END OF CODE ------------------------------
