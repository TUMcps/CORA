function [J_handle,J_str] = jacobianHandle(pZ,id_diff)
% jacobianHandle - computes the function handle for the jacobian of pZ
% with respect to id_diff
%
% Syntax:  
%    PZ = jacobianHandle(pZ)
%    PZ = jacobianHandle(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id - id with respect to which pZ is differentiated
%
% Outputs:
%    J_handle - function handle with argument x and r (r corresponding to
%               ids of pZ.id not in id) returning the jacobian matrix of pZ
%    J_str    - structure matrix indicating zero entries of J_handle for
%               any x,r
%
% Example: 
%    pZ = polyZonotope(0,eye(2),zeros(2,0),[1,0;1,0;0,1],[1;2;3]);
%    fJ = jacobianHandle(pZ);
%    fJ([-1;-2;0])   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: jacobian

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% ids with respect to which we differentiate
if ~exist('id_diff','var')
    id_diff = pZ.id;
elseif length(pZ.id)<length(id_diff) || ~all(ismember(id_diff,pZ.id))
    error('Some ids not contained in pZ.id');
end
id_param = setdiff(pZ.id,id_diff,'stable');
% remove any constant generators
ind_cut_g = all(pZ.G==0,1) | all(pZ.expMat==0,1);
G = pZ.G(:,~ind_cut_g);
expMat = pZ.expMat(:,~ind_cut_g);
% remove any unnecessary ids
ind_cut_id = all(expMat==0,2);
id = pZ.id(~ind_cut_id);
expMat = expMat(~ind_cut_id,:);
% remove any dimensions that are constant
ind_cut_dim = all(pZ.G==0,2);
n = length(pZ.c);
c = pZ.c(~ind_cut_dim);
G = G(~ind_cut_dim,:);

% find used id_diff
id_d = intersect(id_diff,id,'stable');
% find correspondence id_d <=> id_diff
[~,ii_d,~] = intersect(id_diff,id_d,'stable');
% get used id_param
id_p = setdiff(id,id_d,'stable');
% correspondence id_param <=> id_p
[~,ii_p,~] = intersect(id_param,id_p,'stable');
pZ_r = polyZonotope(c,G,pZ.Grest,expMat,id);
%% compute Jacobian
pZ_diff_cell = jacobian(pZ_r,id_d);
d_r = length(id_d);
n_r = length(c);
J_cell = cell(d_r,1);
J_structure = ones(n_r,d_r);
for i=1:d_r
    J_cell{i} = fhandle(pZ_diff_cell{i},{id_d,id_p});
    J_str_i = ~isZero(pZ_diff_cell{i});
    assert(any(J_str_i));
    J_structure(:,i) = J_str_i;
end
%% re-introduce missing dimensions, ids
d = length(id_diff);
J_str = zeros(n,d);
% reintroduce missing dims (additional row vectors), ids (additional
% columns) into J_str
J_str(~ind_cut_dim,ii_d) = J_structure;
% reintroduce missing dims, ids into handle
horzcatc = @(x) horzcat(x{:});
% pad rows to reintroduce missing dimensions; pad cols to reintroduce
% missing ids
J_handle = @(x,p) padZeros(horzcatc(arrayfun(...
                  @(ii) J_cell{ii}(x(ii_d),p(ii_p)),1:d_r,'Uni',false)),...
                  ~ind_cut_dim,ii_d);
if isempty(id_param)
    J_handle = @(x) J_handle(x,[]);
end
end
%%---%%
function X_pad = padZeros(X,rows,cols)
n = length(rows);
m = length(cols);
X_pad = zeros(n,m);
if isa(X,'sym')
    X_pad = sym(X_pad);
end
X_pad(rows,cols) = X;
end
%------------- END CODE --------------