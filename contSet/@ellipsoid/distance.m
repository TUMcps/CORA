function val = distance(E,S)
% distance - computes the distance between an ellipsoid and the another set
%    representation or a point
%
% Syntax:
%    val = distance(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object, numeric, cell-array
%
% Outputs:
%    val - distance(s) between ellipsoid and set/point
%
% Example:
%    E = ellipsoid(eye(2));
%    P = polytope([1 1]/sqrt(2),-2);
%    distance(E,P)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   18-March-2021 (allowing cell arrays)
%                04-July-2022 (VG, replace cell arrays by class arrays)
%                05-October-2024 (MW, remove class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[E,S] = reorderNumeric(E,S);

% check input arguments
inputArgsCheck({{E,'att','ellipsoid'};
                {S,'att',{'contSet','numeric','cell'}}});

% check equal dimensions
equalDimCheck(E,S);

if isnumeric(S)
    val = priv_distancePoint(E,S);
    return
end

% rewrite S as cell-array for easier handling
if ~iscell(S)
    S = {S};
end

% loop over individual pairs
val = zeros(1,numel(S));
for i=1:numel(S)
    val(i) = aux_distance(E,S{i});
end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_distance(E,S)

if representsa_(S,'emptySet',eps)
    % distance to empty set = 0 since empty-set \subseteq obj
    val = 0;
    return;
end

% different distances
if isa(S,'ellipsoid')
    val = priv_distanceEllipsoid(E,S);
    return
end
    
if isa(S,'polytope')
    if representsa_(S,'hyperplane',eps)
        val = priv_distanceHyperplane(E,S);
    else
        val = priv_distancePolytope(E,S);
    end
    return
end

throw(CORAerror('CORA:noops',E,S));

end

% ------------------------------ END OF CODE ------------------------------
