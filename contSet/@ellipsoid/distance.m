function val = distance(E,S)
% distance - Computes the distance between an ellipsoid and the another set
%    representation or a point
%
% Syntax:
%    val = distance(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object/array
%
% Outputs:
%    val - distance(s) between ellipsoid and set/point
%
% Example:
%    E = ellipsoid(eye(2));
%    S = halfspace([1 1]/sqrt(2),-2);
%    distance(E,S)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   18-March-2021 (allowing cell arrays)
%                04-July-2022 (VG, replace cell arrays by class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% parsing and checking
% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}});

% check equal dimensions
equalDimCheck(E,S);

if isa(S,'contSet')
    val = zeros(size(S));
else
    val = zeros(1,size(S,2));
end
if (isnumeric(S) && isempty(S)) || (~isnumeric(S) && all(representsa_(S,'emptySet',eps)))
    % distance to empty set = 0 since empty-set \subseteq obj
    val(:) = 0;
    return;
end

%% different distances
if isa(S,'ellipsoid')
    for i=1:numel(S)
        val(i) = distanceEllipsoid(E,S(i));
    end
    
elseif isa(S,'conHyperplane')
    % conHyperplane actually is a hyperplane
    for i=1:numel(S)
        if representsa_(S(i),'hyperplane',eps)
            val(i) = distanceHyperplane(E,S(i));
        else
            % use polytope implementation
            val(i) = distancePolytope(E,polytope(S(i)));
        end
    end
    
elseif isa(S,'halfspace')
    % convert to polytope
    for i=1:numel(S)
        val(i) = distancePolytope(E,polytope(S(i)));
    end
    
elseif isa(S,'polytope')
    for i=1:numel(S)
        val(i) = distancePolytope(E,S(i));
    end
    
elseif isa(S,'double')
    val = distancePoint(E,S);

else
    throw(CORAerror('CORA:noops',E,S));
    
end

% ------------------------------ END OF CODE ------------------------------
