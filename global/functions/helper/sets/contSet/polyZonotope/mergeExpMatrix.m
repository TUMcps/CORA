function [id,E1,E2] = mergeExpMatrix(id1,id2,E1,E2)
% mergeExpMatrix - Merge the ID-vectors of two polyZonotope objects
%                  and adapte the exponent matrices accordingly
%
% Syntax:
%    [id,E1,E2] = mergeExpMatrix(id1,id2,E1,E2)
%
% Inputs:
%    id1 - ID-vector of the first polynomial zonotope
%    id2 - ID-vector of the second polynomial zonotope
%    E1 - exponent matrix of the first polynomial zonotope
%    E2 - exponent matrix of the second polynomial zonotope
%
% Outputs:
%    id - merged ID-vector
%    E1 - adapted exponent matrix of the first polynomial zonotope
%    E2 - adapted exponent matrix of the second polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       25-June-2018 
% Last update:   07-October-2025 (TL, optimization)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure uniqueness
[E1,id1] = removeRedundantIds(E1,id1);
[E2,id2] = removeRedundantIds(E2,id2);

% gather all ids
id = unique([id1;id2]);

% extend exponent matrices
E1 = (id1' == id) * E1;
E2 = (id2' == id) * E2;

end

% ------------------------------ END OF CODE ------------------------------
