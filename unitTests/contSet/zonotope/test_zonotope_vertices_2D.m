function res = test_zonotope_vertices_2D
% test_zonotope_vertices_2D - unit test function of 2D zonotope/vertices
%    (this was formerly known as zonotope/polygon)
%
% Syntax:
%    res = test_zonotope_vertices_2D
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   11-October-2024 (TL, renamed vertices_2D)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty zonotope
Z = zonotope.empty(2);
p = vertices(Z);
assert(isempty(p) && isnumeric(p) && all(size(p) == [2,0]));

for i=1:2
    if i == 1
        % non-degenerate zonotope
        c = [1;-1];
        G = [2 -3 1 0; -1 1 0 2];
        p_true = [-5 3; -3 3; 3 1; 7 -1; 7 -5; 5 -5; -1 -3; -5 -1]';
    elseif i == 2
        % degenerate zonotope
        c = [2;-1];
        G = [1 0; 0 0];
        p_true = [1 -1; 3 -1]';
    end
    Z = zonotope(c,G);

    % compute polygon points
    p = vertices(Z);
    
    % compare to true points
    
    % true points are a subset because polygon is closed
    assert(compareMatrices(p_true,p,0,'subset'));
    % no differing points
    assert(compareMatrices(p_true,unique(p','rows')'));
    
    % check whether ordering is correct
    angles = atan2d(p(2,:),p(1,:));
    [~,idxMax] = max(angles);
    angles_ = [angles(idxMax+1:end), angles(1:idxMax)];
    assert(all(diff(angles_) >= 0));

end

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
