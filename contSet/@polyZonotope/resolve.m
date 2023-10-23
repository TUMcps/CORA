function res = resolve(pZ,x,varargin)
% resolve - computes result of inserting a value for the identifiers into a
%    polynomial zonotope
%
% Syntax:
%    res = resolve(pZ,x)
%    res = resolve(pZ,x,id)
%
% Inputs:
%    pZ  - polyZonotope object to be resolved
%    x   - value substituted
%    id  - (optional) ids corresponding to x
%
% Outputs:
%    res - resulting polynomial zonotope or double if id==pZ.id
%
% Example: 
%    pZ = polyZonotope([1;-1],eye(2),zeros(2,0),[1,2;0,3]);
%    
%    resolve(pZ,[0.5;1]) % ans = [1.5;-0.75]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: subs

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default value
id = pZ.id;

% parse input arguments
if nargin == 3 && ~isempty(varargin{1})
    if ~isnumeric(varargin{1}) || ~isvector(varargin{1}) || ...
            any(mod(varargin{1},1) ~= 0)
        throw(CORAerror('CORA:wrongValue','third',...
            'Identifiers from polynomial zonotope'));
    end
    id = varargin{1};
end

% check input arguments
if length(id)~=length(x)
    throw(CORAerror('CORA:wrongValue','third',...
        'length of id should match length of input'));
end
if sum(ismember(pZ.id,id))~=length(id) || length(unique(id))~=length(id)
    throw(CORAerror('CORA:wrongValue','third','Invalid id'));
end


%% sort in ascending order such that ismember later on is enough
[id,indid] = sort(id);
x = x(indid);
[Id,ind] = sort(pZ.id);
eM = pZ.E(ind,:);
%% resolve pZ
indx = ismember(Id,id);
eMx = eM(indx,:);
G = pZ.G;
if ~isempty(eMx)
    G = pZ.G.*prod(x.^eMx,1);
end
E = eM;
E(indx,:) = [];
[E,G,c] = removeZeroExponents(E,G);
c = c + pZ.c;
if isempty(E) || all(E(:)==0)
    c = c + sum(G,2);
    G = [];
    E = [];
    if isempty(pZ.GI) || all(pZ.GI(:)==0)
        res = c;
        return;
    end
end
[E,G] = removeRedundantExponents(E,G);

%% revert sort order
ids = pZ.id(~ismember(pZ.id,id));
[~,idx] = sort(ids);
indr_rev(idx) = (1:length(idx));
res = polyZonotope(c, G, pZ.GI, E(indr_rev,:), ids);

% ------------------------------ END OF CODE ------------------------------
