function pZ = restoreId(pZ,id)
% restoreId - adds additional identifiers to a polynomial zonotope; skip
%    identifiers that are already contained in the identifer list
%
% Syntax:
%    res = restoreId(pZ,id)
%
% Inputs:
%    pZ - polyZonotope object
%    id - ids to be restored
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:
%    pZ = polyZonotope(1,1,[],1,2);
%    pZ_ext = restoreId(pZ,[3;2;10]);
%    pZ_ext.id
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: replaceId

% Authors:       Victor Gassmann
% Written:       13-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% exclude identifiers that are already part of the identifier list
ind = ismember(id,pZ.id);
n_toRestore = sum(~ind);

% number of columns in the exponent matrix
ne = size(pZ.E,2);

% append zero-rows to the exponent matrix
pZ.E = [pZ.E;zeros(n_toRestore,ne)];

% append identifiers
pZ.id = [pZ.id;id(~ind)];

% ------------------------------ END OF CODE ------------------------------
