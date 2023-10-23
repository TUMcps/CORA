function [H_handle,H_str] = hessianHandle(pZ,varargin)
% hessianHandle - computes the function handle for the hessian of a
%    one-dimensional polynomial zonotope with respect to id
%
% Syntax:
%    [H_handle,H_str] = hessianHandle(pZ)
%    [H_handle,H_str] = hessianHandle(pZ,tol)
%    [H_handle,H_str] = hessianHandle(pZ,id_diff,id_param)
%    [H_handle,H_str] = hessianHandle(pZ,id_diff,id_param,tol)
%
% Inputs:
%    pZ - polyZonotope object (one-dimensional only)
%    id_diff  - id with respect to which pZ is differentiated
%    id_param - id which should not be differentiated
%    tol-       tolerance
%
% Outputs:
%    H_handle - function handle with argument x (id_diff) and r (id_param) 
%               returning the hessian matrix of pZ
%    H_str    - structure matrix indicating zero entries of H_handle for
%               any x,r
%
% Example: 
%    pZ = polyZonotope(0,[1,2],[],[2,0;1,0;0,1],[1;2;3]);
%    % visualize pZ
%    print(pZ,'Ids',{[1;3],2},'Vars',{'x','p'});
%    fH = hessianHandle(pZ,[1;3],2);
%    % use symbolic toolbox to visualize result
%    x = sym('x',[2,1],'real');
%    p = sym('p',[1,1],'real');
%    fH(x,p)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: jacobianHandle

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   29-November-2021
%                25-February-2022
%                04-July-2022 (VG, moved id-check to common function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(varargin)
    tol = 0;
    [id_diff,id_param] = checkDiffParamIds(pZ);
elseif length(varargin)==1
    tol = varargin{1};
    [id_diff,id_param] = checkDiffParamIds(pZ);
elseif length(varargin)==2
    [id_diff,id_param] = checkDiffParamIds(pZ,varargin{:});
    tol = 0;
elseif length(varargin)==3
    [id_diff,id_param] = checkDiffParamIds(pZ,varargin{1:2});
    tol = varargin{3};
else
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

inputArgsCheck({{tol,'att',{'double'},{'scalar','nonnegative','nonnan'}}});

if dim(pZ) ~= 1
    throw(CORAerror('CORA:wrongValue','first','be one-dimensional'));
end

d = length(id_diff);
ind_d = ismember(pZ.id,id_diff);

% remove any constant or "linear" generators
ind_cut_g = all(withinTol(pZ.G,0,tol),1) | sum(pZ.E(ind_d,:),1)<=1;
G = pZ.G(:,~ind_cut_g);
E = pZ.E(:,~ind_cut_g);

% remove any unnecessary ids
ind_cut_id = all(E==0,2);
id = pZ.id(~ind_cut_id);
E = E(~ind_cut_id,:);

% find used id_diff
id_d = id_diff(ismember(id_diff,id));

pZ.G = G; pZ.E = E; pZ.id = id;
if isempty(pZ.G)
    H_str = zeros(d);
    if isempty(id_param)
        H_handle = @(x) zeros(d);
    else
        H_handle = @(x,p) zeros(d);
    end
    return;
end

%% first and second derivative
% id_d in occurence of id_diff
[id_diff_d,ii_diff_d] = intersect(id_diff,id_d,'stable');
pZ_diff_cell = jacobian(pZ,id_diff_d);
H_cell = repmat({@(x,p)0},d,d);
H_str = zeros(d);

for i=1:length(pZ_diff_cell)
    index_col = ii_diff_d(i);
    pZd_i = pZ_diff_cell{i};
    % remove redundant generators (constant)
    ind_cut_gi = withinTol(pZd_i.G,0,tol) | sum(pZd_i.E,1)==0;
    Gi = pZd_i.G(:,~ind_cut_gi);
    eMi = pZd_i.E(:,~ind_cut_gi);
    % remove redundant ids
    ind_cut_i = all(eMi==0,2);
    % used ids
    id_i = pZd_i.id(~ind_cut_i);
    eMi = eMi(~ind_cut_i,:);
    % used diff ids
    id_di = id_diff(ismember(id_diff,id_i));
    assert(~isempty(id_di),'Further reduction possible!');
    pZ_diff_i = polyZonotope(0*pZd_i.c,Gi,pZd_i.GI,eMi,id_i); % 0*.. test
    assert(~representsa_(pZ_diff_i,'origin',tol),'Bug: Can be simplified further beforehand!');
    
    % get indices at which to start computation (due to symmetry of hessian)
    % get id_di as it occurs in id_diff 
    [~,ii_diff_di,~] = intersect(id_diff,id_di,'stable');
    ind_rem = ii_diff_di>=i;
    if ~any(ind_rem)
        continue;
    end
    ii_diff_rem = ii_diff_di(ind_rem);
    id_diff_rem = id_diff(ii_diff_rem);
    pZHi_cell = jacobian(pZ_diff_i,id_diff_rem); 
    for j=1:length(id_diff_rem)
        index_row = ii_diff_rem(j);
        %assert(~representsa_(pZHi_cell{j},'origin',eps),'Bug: Should be true?');
        H_str(index_row,index_col) = 1;
        % find which ids remain of id_diff and id_param
        id_dij = id_diff(ismember(id_diff,pZHi_cell{j}.id));
        id_pij = id_param(ismember(id_param,pZHi_cell{j}.id));
        assert(all(ismember([id_dij;id_pij],pZHi_cell{j}.id)),'Bug: should contain all ids');
        % find correspondence
        ii_dij = subsetIndex(id_diff,id_dij);
        if isempty(id_pij)
            ftmp = fhandle(pZHi_cell{j},{id_dij});
            % reintroduce previously removed ids
            H_cell{index_row,index_col} = @(x,p) ftmp(x(ii_dij));
        else
            ii_pij = subsetIndex(id_param,id_pij);
            ftmp = fhandle(pZHi_cell{j},{id_dij,id_pij});
            % reintroduce previously removed ids
            H_cell{index_row,index_col} = @(x,p) ftmp(x(ii_dij),p(ii_pij));
        end
        
    end
end

%% construct handle
symmetric = @(M) tril(M) + triu(M',1);
H_handle = @(x,p) symmetric(aux_symwrap(cellfun(@(f)f(x,p),H_cell,'Uni',false)));
H_str = symmetric(H_str);
if isempty(id_param)
    H_handle = @(x) H_handle(x,[]);
end

end


% Auxiliary functions -----------------------------------------------------

function M = aux_symwrap(M)

    if any(cellfun(@(cc)isa(cc,'sym'),M(:)))
        M = cellfun(@(cc)sym(cc),M);
    else
        M = cellfun(@(cc)cc,M);
    end

end

% ------------------------------ END OF CODE ------------------------------
