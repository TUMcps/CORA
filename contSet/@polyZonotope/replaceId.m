function pZ = replaceId(pZ,varargin)
% replaceId - replace all or some of the identifiers with new ones
%
% Syntax:
%    pZ = replaceId(pZ,id_part_new)
%    pZ = replaceId(pZ,id_part_old,id_part_new)
%
% Inputs:
%    pZ - polyZonotope object
%    id_part_new - new identifiers replacing old identifiers
%    id_part_old - old identifers
%
% Outputs:
%    pZ - modified polyZonotope object with replaced identifiers
%
% Example:
%    pZ = polyZonotope(0,[1,2],[],[1,2;3,4;5,6],[1;2;3]);
%    pZ_newId = replaceId(pZ,[3;1],[1;3]);
%    pZ_newId.id
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       12-January-2021 
% Last update:   28-February-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if length(varargin) == 1
    id_part_old = pZ.id;
    id_part_new = varargin{1};
elseif length(varargin) == 2
    id_part_old = varargin{1};
    id_part_new = varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% check input arguments
if length(id_part_old)~=length(id_part_new) || ~all(ismember(id_part_old,pZ.id))
     throw(CORAerror('CORA:wrongValue','second and/or third',...
         "The identifiers in 'id_part_new' and 'id_part_old' should be of the same length, " + ...
         "and they should only contain identifiers contained in the polynomial zonotope."));
end

% read out indices of identifiers in polynomial zonotope
ii_o = subsetIndex(pZ.id,id_part_old);
pZ.id(ii_o) = id_part_new;

% remove possibly duplicate entries from id
[id_new,~,iic] = unique(pZ.id,'stable'); 
ii_id = accumarray(iic,(1:length(pZ.id))',[],@(x){x});
E = zeros(length(ii_id),size(pZ.E,2));
for i=1:length(ii_id)
    E(i,:) = sum(pZ.E(ii_id{i},:),1);
end

% assign identifiers and exponent matrix
pZ.id = id_new;
pZ.E = E;

% ------------------------------ END OF CODE ------------------------------
