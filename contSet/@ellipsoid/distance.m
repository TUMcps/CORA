function Val = distance(obj,S)
% enclose - Computes the distance between an ellipsoid and the set S (or
% single point)
%
% Syntax:  
%    E = distance(obj,S)
%
% Inputs:
%    obj - ellipsoid object
%    S   - set (ellipsoid, constHyperplane, ...) or cell array thereof
%
% Outputs:
%    Val - distance(s) between obj and S
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  18-March-2021 (allowing cell arrays)
% Last revision:---

%------------- BEGIN CODE --------------
%% parsing & checking
if ~isa(obj,'ellipsoid')
    error('First input argument has to be of type "ellipsoid"!');
end

% handle empty cells, S not a cell etc
[S,size_S] = prepareSetCellArray(S,obj);
if isempty(S)
    % distance to empty set = 0 since empty-set \subseteq obj
    Val = 0;
    return;
end
N = length(S);
%% different distances
Val = zeros(N,1);
if isa(S{1},'ellipsoid')
    for i=1:N
        Val(i) = distanceEllipsoid(obj,S{i});
    end
    
elseif isa(S{1},'conHyperplane')
    % conHyperplane actually is a hyperplane
    for i=1:N
        if isHyperplane(S{i})
            Val(i) = distanceHyperplane(obj,S{i});
        else
            % use mptPolytope implementation
            S{i} = mptPolytope(S{i});
            Val(i) = distanceMptPolytope(E,S{i});
        end
    end
    
elseif isa(S{1},'halfspace')
    % convert to mptPolytope
    for i=1:N
        Val(i) = distanceMptPolytope(obj,mptPolytope(S{i}));
    end
    
elseif isa(S{1},'mptPolytope')
    for i=1:N
        Val(i) = distanceMptPolytope(obj,S{i});
    end
    
elseif isa(S{1},'double')
    Val = distancePoint(obj,S);
else
    error('second input argument type not supported');
end
Val = reshape(Val,size_S);
%------------- END OF CODE --------------