function [H,H_str] = addedTermsHessian(pZ,varargin)
% addedTermsHessian - computes a weighted addition using all dimensions of 
%                     polyZonotope and then calculates the hessian handle
%                     for the result
%
%
% Syntax:
%    [H,H_str] = addedTermsHessian(pZ)
%    [H,H_str] = addedTermsHessian(pZ,tol)
%    [H,H_str] = addedTermsHessian(pZ,id_diff,id_param)
%    [H,H_str] = addedTermsHessian(pZ,id_diff,id_param,tol)
%
% Inputs:
%    pZ      - polyZonotope object
%    id_diff - ids of pZ.id wrt which we compute hessian
%    id_param- all remaining ids in pZ.id (assumed to be fixed)
%    tol     - tolerance
%
% Outputs:
%    H      - hessian handle: H(x,p,l) or H(x,l) (if id_param = []) where
%                        H(x,p,l) = @(x,p,l)
%                        hessian_x(l(1)*pZ_1(x,p)+l(2)*pZ_2(x,p)+...)
%    H_str  - structure matrix indicating zero entries of H_handle for any
%             x,p,l
%
% Example: 
%    pZ = polyZonotope([1;2;3],[2,1,1,0,0;0,0,0,1,0;0,0,0,0,1],[],...
%                      [0,2,0,1,0;0,1,1,1,2;1,0,1,0,0;0,1,1,0,2],...
%                       [5;2;-1;-2]);
%    % visualize pZ (idx_1 = 5, idx_2 = 2, idp_1 = -1; idp_2 = -2;
%    print(pZ,'Ids',{[5;2],[-1;-2]},'Vars',{'x','p'});
%    % compute weighted sum hessian
%    [f_H,H_str] = addedTermsHessian(pZ,[5;2],[-1;-2]);
%    % use symbolic toolbox to check result
%    x = sym('x',[2,1],'real');
%    p = sym('p',[2,1],'real');
%    l = sym('l',[max(abs(pZ.id))+dim(pZ),1],'real');
%    f_H(x,p,l)
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       22-December-2021
% Last update:   04-July-2022 (VG, error handeling, additional parameters)
%                16-November-2022 (LS, return structure matrix)
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

% number of terms to add and weight
Nt = dim(pZ);

id = pZ.id;

ind_diff = ismember(id,id_diff);

% cut all terms that would vanish if differentiated twice
ind_cut = sum(pZ.E(ind_diff,:),1)<=1 | all(withinTol(pZ.G,0,tol),1);
if ~all(ind_cut)
    pZ_gsum = [];
    % ids for lambda
    idl = (max(abs(id))+(1:Nt))';
    for i=1:Nt
        pZ_i = project(pZ,i);
        % again, cut all terms that would vanish if differentiated twice
        ind_cut_i = withinTol(pZ_i.G,0,tol) | sum(pZ_i.E(ind_diff,:),1)<=1;

        % remove <= linear generators
        if ~all(ind_cut_i)
            pZ_i.G = pZ_i.G(:,~ind_cut_i);
            pZ_i.E = pZ_i.E(:,~ind_cut_i);

            pZ_i_ext = pZ_i;
            % add dep factor used to input weight later on (will not be
            % differentiated against)
            mi = size(pZ_i.E,2);
            pZ_i_ext.E = [pZ_i.E;ones(1,mi)];
            pZ_i_ext.id = [pZ_i.id;idl(i)];
            % remove center (since it will be thrown away anyway)
            pZ_i_ext.c = 0;
            % sum to result
            if ~isempty(pZ_gsum)
                pZ_gsum = exactPlus(pZ_gsum,pZ_i_ext);
            else
                pZ_gsum = pZ_i_ext;
            end
        end
    end
    
    % compute which ids "survived" (also necessary below to recover
    % "surviving" ids for function handle)
    ind_diff_rem = ismember(id_diff,pZ_gsum.id);
    ind_param_rem = ismember(id_param,pZ_gsum.id);
    ind_l_rem = ismember(idl,pZ_gsum.id);
    id_diff_rem = id_diff(ind_diff_rem);
    idl_rem = idl(ind_l_rem);
    id_param_rem = id_param(ind_param_rem);
    % compute hessian of summed average where we only differentiate with
    % respect to id_diff_rem, and treat id_param_rem and idl_rem as fixed
    % parameters
    [H_ctr,H_str] = hessianHandle(pZ_gsum,id_diff_rem,[id_param_rem;idl_rem],tol);
    % reintroduce the "lost" x's, p's, l's
    if isempty(id_param)
        H = @(x,l) H_ctr(x(ind_diff_rem),l(ind_l_rem));
    else
        H = @(x,p,l) H_ctr(x(ind_diff_rem),[p(ind_param_rem);l(ind_l_rem)]);
    end

else
    % linear (or constant) constraints
    if isempty(id_param)
        H = @(x,l) zeros(length(id_diff));
    else
        H = @(x,p,l) zeros(length(id_diff));
    end
    H_str = zeros(length(id_diff));
end

% ------------------------------ END OF CODE ------------------------------
