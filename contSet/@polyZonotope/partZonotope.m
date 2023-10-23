function [PZ,pZ_gi_cell,ind_constGen,EM,indexE] = partZonotope(pZ,id)
% partZonotope - computes a zonotope overapproximation in id
%
% Syntax:
%    [PZ,pZ_gi_cell] = partZonotope(pZ,id)
%    [PZ,pZ_gi_cell,ind_constGen,EM,indexE] = partZonotope(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id - ???
%
% Outputs:
%    PZ - zonotope overapproximation of pZ in id (in general a polyZonotope
%                 object)
%    pZ_gi_cell - cell array of polyZonotopes (generators of new zonotope
%                 factors)
%    ind_constGen - ???
%    EM - ???
%    indexE - ???
%
% Example: 
%    pZ = polyZonotope(1,[1,1],[],[2,0;2,1;0,1],[1;2;3]);
%    partZono = partZonotope(pZ,2);
%    id_new = partZono.id(~ismember(partZono.id,pZ.id));
%    print(partZono,'Ids',{[1;3],id_new},'Vars',{'b','w'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Victor Gassmann
% Written:       24-March-2020
% Last update:   07-March-2022 (simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% overapproximates pZ as a zonotope only in id_part
% pZ_gi_cell: generators of {pZ}_w which are non-constant in p
if length(id)>length(pZ.id) || ~all(ismember(id,pZ.id))
    throw(CORAerror('CORA:wrongValue','second',...
        "The length of 'id' should not exceed that of pZ.id, " + ...
        "and id should only contain identifiers contained in pZ.id."));
end

nx = dim(pZ);
idv = id;
idp = pZ.id(~ismember(pZ.id,idv));
G = [pZ.c,pZ.G];
E = [zeros(size(pZ.E,1),1),pZ.E];
ind_v = ismember(pZ.id,idv);
ind_p = ismember(pZ.id,idp);
np = length(idp);

%% compute overapproximation
eM_v = E(ind_v,:);
[~,~,Indc] = unique(eM_v','rows','stable');
indexE = accumarray(Indc,(1:size(E,2))',[],@(x){x});
% first is all-zero vector by design (see l.45, l.47)
% extract all constant-in-v generators (includes center of course)
G_0 = G(:,indexE{1});
eMp_0 = E(ind_p,indexE{1});
[eMp_0,G_0,c_0] = removeZeroExponents(eMp_0,G_0);
% center zonotope only dependent on constants and p
pZ_cp = polyZonotope(c_0,G_0,zeros(nx,0),eMp_0,idp);

pZ_gi_cell = cell(length(indexE),1);
ind_constGen = false(1,length(indexE));
EM = zeros(length(idv),length(indexE));
% first id is not used since start at 2 (center already handled)
ido = (max([idp;idv])+1:max([idp;idv])+length(indexE)-1)';

PZ = polyZonotope(zeros(nx,1),zeros(nx,0),zeros(nx,0),zeros(length(ido),0),ido);
for i=2:length(indexE)
    % unique v exponent vector
    e_vi = E(ind_v,indexE{i}(1));
    EM(:,i) = e_vi;
    % all exponents with same v part
    eMp_i = E(ind_p,indexE{i});
    Gi = G(:,indexE{i});
    [eMp_i,Gi,ci] = removeZeroExponents(eMp_i,Gi);
    if sum(mod(e_vi,2))==0
        % only even exponents
        pZtmp = polyZonotope(1/2*ci,1/2*Gi,zeros(nx,0),eMp_i,idp);
        pZ_cp = exactPlus(pZ_cp,pZtmp);
        ci = 1/2*ci;
        Gi = 1/2*Gi;
    end
    
    tmp = polyZonotope(zeros(nx,1),[ci,Gi],zeros(nx,0),...
                [zeros(np,1),eMp_i;ones(1,1+size(eMp_i,2))],[idp;ido(i-1)]);
    % add generator (center added at the end)
    PZ = exactPlus(PZ,tmp);
    % check if e_vi corresponds only to constant-in-p generators
    if isempty(eMp_i)
        ind_constGen(i) = true;
        pZ_gi_cell{i} = polyZonotope(ci,zeros(nx,0),zeros(nx,0),zeros(np,0),idp);
    else
        pZ_gi = polyZonotope(ci,Gi,zeros(nx,0),eMp_i,idp);
        pZ_gi_cell{i} = pZ_gi; 
    end
end
% add center
PZ = exactPlus(PZ,pZ_cp);
% add rest matrix
PZ.GI = pZ.GI;

% order important (center first important!)
pZ_gi_cell{1} = pZ_cp;

% ------------------------------ END OF CODE ------------------------------
