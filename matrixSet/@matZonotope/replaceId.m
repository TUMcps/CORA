function mZ = replaceId(mZ,varargin)
% replaceId - replace all or some of the identifiers with new ones
%
% Syntax:
%    pZ = replaceId(mZ,id_part_new)
%    pZ = replaceId(mZ,id_part_old,id_part_new)
%
% Inputs:
%    pZ - matZonotope object
%    id_part_new - new identifiers replacing old identifiers
%    id_part_old - old identifers
%
% Outputs:
%    mZ - modified matZonotope object with replaced identifiers
%
% Example:
%    mZ = polyZonotope(0,[1,2],[],[1,2;3,4;5,6],[1;2;3]);
%    mZ_newId = replaceId(mZ,[3;1],[1;3]);
%    mZ_newId.id
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Laura Luetzow
% Written:       25-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,3);

% parse input arguments
if length(varargin) == 1
    id_part_old = mZ.id;
    id_part_new = varargin{1};
elseif length(varargin) == 2
    id_part_old = varargin{1};
    id_part_new = varargin{2};
end

% check input arguments
if length(id_part_old)~=length(id_part_new) || ~all(ismember(id_part_old,mZ.id))
     throw(CORAerror('CORA:wrongValue','second and/or third',...
         "The identifiers in 'id_part_new' and 'id_part_old' should be of the same length, " + ...
         "and they should only contain identifiers contained in the matrix zonotope."));
end

% read out indices of identifiers in matrix zonotope
ii_o = subsetIndex(mZ.id,id_part_old);
mZ.id(ii_o) = id_part_new;

% ------------------------------ END OF CODE ------------------------------
