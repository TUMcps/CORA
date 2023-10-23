function [pZ_id,pZ_r] = onlyId(pZ,id)
% onlyId - splits pZ into part only dependent on 'id' and rest
%
% Syntax:
%    [pZ_id,pZ_r] = onlyId(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id - identifiers
%
% Outputs:
%    pZ_id - polyZonotope containing only the center and generators of pZ
%            which only contain ids from "id"
%    pZ_r  - all generators from pZ_id not included in pZ_id (and rest
%            matrix)
%
% Example: 
%    pZ = polyZonotope([1;-2],[1,2,3;4,5,6],0.1*eye(2),[1,1,0;1,0,1],[1;2]);
%    print(pZ,'Ids',{1,2},'Vars',{'x','y'});
%    [pZ_x,pZ_xy] = onlyId(pZ,1);
%    print(pZ_x,'Ids',{1},'Vars',{'x'});
%    print(pZ_xy,'Ids',{1,2},'Vars',{'x','y'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{id,'att','double',{'integer','ncols',1}}});

% check if unique
if length(unique(id))<length(id)
    throw(CORAerror('CORA:wrongValue','second',...
        "'id' should not contain any duplicates."))
end

% check if 'id' contains only ids from pZ.id
if ~all(ismember(id,pZ.id))
    throw(CORAerror('CORA:wrongValue','second',...
        "'id' should only contain identifers of the polynomial zonotope."));
end


ind_id = ismember(pZ.id,id);
n = dim(pZ);
% extract all exponent vectors which only have ids from 'id'
ind_mask = sum(pZ.E(~ind_id,:),1)==0;
G_id = pZ.G(:,ind_mask);
eM_id = pZ.E(ind_id,ind_mask);
% put everything containing also ids not in 'id' into rest polyZonotope
G_r = pZ.G(:,~ind_mask);
eM_r = pZ.E(:,~ind_mask);
% construct results
pZ_id = polyZonotope(pZ.c,G_id,zeros(n,0),eM_id,pZ.id(ind_id));
pZ_r = polyZonotope(zeros(n,1),G_r,pZ.GI,eM_r,pZ.id);

% ------------------------------ END OF CODE ------------------------------
