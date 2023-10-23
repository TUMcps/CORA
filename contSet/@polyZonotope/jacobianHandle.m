function [J_handle,J_str] = jacobianHandle(pZ,varargin)
% jacobianHandle - computes the function handle for the jacobian of pZ
%    with respect to id_diff
%
% Syntax:
%    PZ = jacobianHandle(pZ)
%    PZ = jacobianHandle(pZ,id_diff,id_param)
%
% Inputs:
%    pZ - polyZonotope object
%    id_diff - id with respect to which pZ is differentiated
%    id_param- not differentiated ids
%
% Outputs:
%    J_handle - function handle with argument x and r (r corresponding to
%               ids of pZ.id not in id) returning the jacobian matrix of pZ
%    J_str - structure matrix indicating zero entries of J_handle for
%               any x,r
%
% Example: 
%    pZ = polyZonotope([1;-2],eye(2),zeros(2,0),[1,0;1,0;0,1],[1;2;3]);
%    % visualize pZ
%    print(pZ,'Ids',{pZ.id},'Vars',{'x'});
%    fJ = jacobianHandle(pZ);
%    % use symbolic toolbox to visualize result
%    x = sym('x',[3,1],'real');
%    fJ(x)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: jacobian

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   29-November-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check id inputs
[id_diff,id_param] = checkDiffParamIds(pZ,varargin{:});

% remove any constant generators
ind = ismember(pZ.id,id_diff);
ind_cut_g = all(pZ.G==0,1) | all(pZ.E(ind,:)==0,1);
G = pZ.G(:,~ind_cut_g);
E = pZ.E(:,~ind_cut_g);

% remove any unnecessary ids
ind_cut_id = all(E==0,2);
id = pZ.id(~ind_cut_id);
E = E(~ind_cut_id,:);

% remove any dimensions that are constant
ind_cut_dim = all(pZ.G==0,2);
n = length(pZ.c);
d = length(id_diff);
c = pZ.c(~ind_cut_dim);
G = G(~ind_cut_dim,:);

% find used id_diff
id_d = intersect(id_diff,id,'stable');
% find correspondence id_d <=> id_diff: id_diff(ii_d) = id_d
ii_d = subsetIndex(id_diff,id_d);

% get used id_param
id_p = intersect(id_param,id,'stable');
% correspondence id_param <=> id_p: id_param(ii_p) = id_p
ii_p = subsetIndex(id_param,id_p);
pZ.c = c; pZ.G = G; pZ.E = E; pZ.id = id;
if isempty(pZ.G)
    J_str = zeros(n,d);
    if isempty(id_param)
        J_handle = @(x) zeros(n,d);
    else
        J_handle = @(x,p) zeros(n,d);
    end
    return;
end

%% compute Jacobian
pZ_diff_cell = jacobian(pZ,id_d);
d_r = length(id_d);
n_r = length(c);
J_cell = cell(d_r,1);
J_structure = ones(n_r,d_r);
for i=1:d_r
    J_cell{i} = fhandle(pZ_diff_cell{i},{id_d,id_p});
    J_str_i = ~representsa_(pZ_diff_cell{i},'origin',0);
    assert(any(J_str_i));
    J_structure(:,i) = J_str_i;
end

%% re-introduce missing dimensions, ids
J_str = zeros(n,d);
% reintroduce missing dims (additional row vectors), ids (additional
% columns) into J_str
J_str(~ind_cut_dim,ii_d) = J_structure;
% reintroduce missing dims, ids into handle
horzcatc = @(x) horzcat(x{:});
% pad rows to reintroduce missing dimensions; pad cols to reintroduce
% missing ids
% convert ii_d to logical mask
ind_id_d = ismember(1:d,ii_d);
J_handle = @(x,p) aux_padZeros(horzcatc(arrayfun(...
                  @(ii) J_cell{ii}(x(ii_d),p(ii_p)),1:d_r,'Uni',false)),...
                  ~ind_cut_dim,ind_id_d);
if isempty(id_param)
    J_handle = @(x) J_handle(x,[]);
end

end


% Auxiliary functions -----------------------------------------------------

function X_pad = aux_padZeros(X,rows,cols)
    
    n = length(rows);
    m = length(cols);
    X_pad = zeros(n,m);
    if isa(X,'sym')
        X_pad = sym(X_pad);
    end
    X_pad(rows,cols) = X;
    
end

% ------------------------------ END OF CODE ------------------------------
