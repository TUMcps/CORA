function res = test_zonotope_polygon
% test_zonotope_polygon - unit test function of polygon
%
% Syntax:
%    res = test_zonotope_polygon
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty zonotope
Z = zonotope.empty(2);
p = polygon(Z);
res = isempty(p) && isnumeric(p) && all(size(p) == [2,0]);

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
    p = polygon(Z);
    
    % compare to true points
    
    % true points are a subset because polygon is closed
    res(end+1,1) = compareMatrices(p_true,p,0,'subset');
    % no differing points
    res(end+1,1) = compareMatrices(p_true,unique(p','rows')');
    % first and last points need to be the same
    res(end+1,1) = all(withinTol(p(:,1),p(:,end)));
    
    % check whether ordering is correct
    angles = atan2d(p(2,:),p(1,:));
    [~,idxMax] = max(angles);
    angles_ = [angles(idxMax+1:end), angles(1:idxMax)];
    res(end+1,1) = all(diff(angles_) >= 0);

end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
