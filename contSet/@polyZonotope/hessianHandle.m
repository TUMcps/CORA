function [H_handle,H_str] = hessianHandle(pZ,id_diff)
% hessianHandle - computes the function handle for the hessian of 1D-pZ
% with respect to id
%
% Syntax:  
%    PZ = hessianHandle(pZ)
%    PZ = hessianHandle(pZ,id_diff)
%
% Inputs:
%    pZ - polyZonotope object
%    id_diff - id with respect to which pZ is differentiated
%
% Outputs:
%    H_handle - function handle with argument x and r (r corresponding to
%               ids of pZ.id not in id) returning the hessian matrix of pZ
%    H_str    - structure matrix indicating zero entries of H_handle for
%               any x,r
%
% Example: 
%    pZ = polyZonotope(0,eye(2),zeros(2,0),[2,0;1,0;0,1],[1;2;3]);
%    fH = hessianHandle(pZ);
%    fH([-1;-2;0])   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: jacobianHandle

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% Note: Fact that hessian is symmetrical is used to reduce
% computational load
% ids with respect to which we differentiate
if ~exist('id_diff','var')
    id_diff = pZ.id;
elseif length(pZ.id)<length(id_diff) || ~all(ismember(id_diff,pZ.id))
    error('Some ids not contained in pZ.id!');
end
n = length(pZ.c);
d = length(id_diff);
if n~=1
    error('pZ has to be of dimension 1!');
end
id_param = setdiff(pZ.id,id_diff,'stable');
% remove any constant or "linear" generators
ind_cut_g = all(pZ.G==0,1) | sum(pZ.expMat,1)<=1;
G = pZ.G(:,~ind_cut_g);
expMat = pZ.expMat(:,~ind_cut_g);
% remove any unnecessary ids
ind_cut_id = all(expMat==0,2);
id = pZ.id(~ind_cut_id);
expMat = expMat(~ind_cut_id,:);

% find used id_diff
id_d = intersect(id_diff,id,'stable');
% find correspondence id_d <=> id_diff
[~,ii_d,~] = intersect(id_diff,id_d,'stable');
pZ_r = polyZonotope(pZ.c,G,pZ.Grest,expMat,id);
%% first & second derivative
pZ_diff_cell = jacobian(pZ_r,id_d);
d_r = length(id_d);
H_cell = repmat({@(x,p)0},d,d);
H_str = zeros(d);
for i=1:d_r
    pZd_i = pZ_diff_cell{i};
    % remove redundant generators (constant)
    ind_cut_gi = pZd_i.G==0 | sum(pZd_i.expMat,1)==0;
    Gi = pZd_i.G(:,~ind_cut_gi);
    eMi = pZd_i.expMat(:,~ind_cut_gi);
    % remove redundant ids
    ind_cut_i = all(eMi==0,2);
    % used ids
    id_i = id(~ind_cut_i);
    eMi = eMi(~ind_cut_i,:);
    % used diff ids
    id_di = intersect(id_i,id_d,'stable');
    pZ_diff_i = polyZonotope(pZd_i.c,Gi,pZd_i.Grest,eMi,id_i);
    assert(~isZero(pZ_diff_i),'Bug: Can be simplified further beforehand!');
    
    % used param ids
    id_pi = setdiff(id_i,id_di,'stable');
    % correspondence between id_diff and id_di
    [~,ii_di,~] = intersect(id_diff,id_di,'stable');
    % correspondence between id_param and id_pi
    [~,ii_pi,~] = intersect(id_param,id_pi,'stable');
    % get indices at which to start computation (due to symmetry of hessian)
    ii_di_rem = ii_di(ii_d(i)<=ii_di);
    pZHi_cell = jacobian(pZ_diff_i,id_diff(ii_di_rem)); 
    for j=1:length(ii_di_rem)
        assert(~isZero(pZHi_cell{j}),'Bug: Should be true?');
        H_str(ii_di_rem(j),ii_d(i)) = 1;
        ftmp = fhandle(pZHi_cell{j},{id_di,id_pi});
        % reintroduce previously removed ids
        H_cell{ii_di_rem(j),ii_d(i)} = @(x,p) ftmp(x(ii_di),p(ii_pi));
    end
end
%% construct handle
symmetric = @(M) tril(M) + triu(M',1);
H_handle = @(x,p) symmetric(symwrap(cellfun(@(f)f(x,p),H_cell,'Uni',false)));
H_str = symmetric(H_str);
if isempty(id_param)
    H_handle = @(x) H_handle(x,[]);
end
end

%---- helper
function M = symwrap(M)
if any(cellfun(@(cc)isa(cc,'sym'),M(:)))
    M = cellfun(@(cc)sym(cc),M);
else
    M = cellfun(@(cc)cc,M);
end
end