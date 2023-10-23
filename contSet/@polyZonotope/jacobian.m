function pZdiff_cell = jacobian(pZ,varargin)
% jacobian - computes the derivatives of each dimension of pZ with respect
%    to id_diff
%
% Syntax:
%    PZ = jacobian(pZ)
%    PZ = jacobian(pZ,id_diff)
%
% Inputs:
%    pZ - polyZonotope object (n-dimensional)
%    id_diff - (optional) ids with respect to which pZ is differentiated
%              (d-dimensional)
%
% Outputs:
%    PZ - d-dimensional cell array where each element is a polyZonotope of
%         dimension n
%
% Example: 
%    pZ = polyZonotope([1;2],[1,2;3,4],[],[1,1;3,2;0,1],[1;2;3]);
%    pZ_jac_cell = jacobian(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cartProd

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default value
id_diff = pZ.id;

% parse input arguments
if nargin == 2 && ~isempty(varargin{1})
    if ~isnumeric(varargin{1}) || ~isvector(varargin{1}) || ...
            any(mod(varargin{1},1) ~= 0) || length(pZ.id)<length(id_diff) || ...
            ~all(ismember(id_diff,pZ.id))
        throw(CORAerror('CORA:wrongValue','second',...
            "All entries in id_diff have to be identifiers in the polynomial zonotope."));
    end
    id_diff = varargin{1};
end

n = dim(pZ);
d = length(id_diff);
I = eye(length(pZ.id));
% will have d-dimensional cell array, since we choose to differentiate the
% n-dimensional polyZonotope for each id to be differentiated
pZdiff_cell = cell(d,1);
% compute correspondence between pZ.id and id_diff: id_diff=pZ.id(ii_diff)
[~,~,ii_diff] = intersect(id_diff,pZ.id,'stable');
for i=1:d
    % for each entry in id_diff, differentiate by subtracting 1 from each
    % column in the appropriate row
    eMtmp = pZ.E - I(:,ii_diff(i));
    % remove any negative entries (happens if constant previously)
    indtmp = any(eMtmp<0,1);
    eMtmp(:,indtmp) = [];
    if ~isempty(eMtmp)
        % adapt G matrix accordingly
        Gtmp = pZ.E(ii_diff(i),~indtmp).*pZ.G(:,~indtmp);
        [E,G,c] = removeZeroExponents(eMtmp,Gtmp);
    else
        E = zeros(length(pZ.id),0);
        G = zeros(n,0);
        c = zeros(n,1);
    end
    % construct result for i-th id_diff entry
    pZdiff_cell{i} = polyZonotope(c,G,zeros(n,0),E,pZ.id);
end

% ------------------------------ END OF CODE ------------------------------
