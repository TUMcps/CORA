function res = resolve(pZ,x,id)
% resolve - computes result of inserting value x for id into pZ
%
% Syntax:  
%    res = resolve(pZ,x)
%    res = resolve(pZ,x,id)
%
% Inputs:
%    pZ  - polyZonotope object to be resolved
%    x   - value substituted
%    id  - ids corresponding to x
%
% Outputs:
%    res - resulting polynomial Zonotope or double if id==pZ.id
%
% Example: 
%    % substitution
%    pZ = polyZonotope([1;-1],eye(2),zeros(2,0),[1,2;0,3]);
%    
%    resolve(pZ,[0.5;1])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: subs

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('id','var')
    id = pZ.id;
end
if length(id)~=length(x)
        error('Length of id (default: all) does not match length of input');
end
if sum(ismember(pZ.id,id))~=length(id) || length(unique(id))~=length(id)
    error('Invalid id');
end
%% sort in ascending order such that ismember later on is enough
[id,indid] = sort(id);
x = x(indid);
[Id,ind] = sort(pZ.id);
eM = pZ.expMat(ind,:);
%% resolve pZ
indx = ismember(Id,id);
eMx = eM(indx,:);
% speed up computation if only identity matrix
if all(size(eMx)==size(eMx,1)) && all(diag(eMx)==1)
    G = pZ.G.*x';
else
    G = pZ.G.*prod(x.^eMx,1);
end
ExpMat = eM;
ExpMat(indx,:) = [];
[ExpMat,G,c] = removeZeroExponents(ExpMat,G);
C = c + pZ.c;
if isempty(ExpMat) || all(ExpMat(:)==0)
    C = C + sum(G,2);
    G = [];
    ExpMat = [];
    if isempty(pZ.Grest) || all(pZ.Grest(:)==0)
        res = C;
        return;
    end
end
[ExpMat,G] = removeRedundantExponents(ExpMat,G);
%% revert sort order
idr = pZ.id(~ismember(pZ.id,id));
[~,tmp1] = sort(idr);
indr_rev(tmp1) = (1:length(tmp1));
res = polyZonotope(C,G,pZ.Grest,ExpMat(indr_rev,:),idr);
%------------- END OF CODE --------------