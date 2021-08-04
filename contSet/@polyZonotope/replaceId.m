function PZ = replaceId(pZ,varargin)
% replaceId - replace (part of) id with new one
%
% Syntax:  
%    PZ = replaceId(pZ,id_new)
%    PZ = replaceId(pZ,id_old,id_new)
%
% Inputs:
%    pZ     - polyZonotope object
%    id_new - new id that replaces id_old
%    id_old - old id
%
% Outputs:
%    PZ     - polyZonotope object with id_old replaced by id_new
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Victor Gassmann
% Written:       12-January-2021 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
if length(varargin)==1
    id_part_old = pZ.id;
    id_part_new = varargin{1};
elseif length(varargin)==2
    id_part_old = varargin{1};
    id_part_new = varargin{2};
else
    error('Wrong type of input arguments');
end
if length(id_part_old)~=length(id_part_new) || ~all(ismember(id_part_old,pZ.id))
    error('id not valid');
end
[Id,ind] = sort(pZ.id);
[Id_part_old,ind_part_old] = sort(id_part_old);
Id_part_new = id_part_new(ind_part_old);
Id_new = Id;
Id_new(ismember(Id,Id_part_old)) = Id_part_new;
ind_rev(ind) = 1:length(ind);
PZ = polyZonotope(pZ.c,pZ.G,pZ.Grest,pZ.expMat,Id_new(ind_rev));
%------------- END OF CODE --------------