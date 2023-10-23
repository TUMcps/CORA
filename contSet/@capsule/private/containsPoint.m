function res = containsPoint(C,p)
% containsPoint - determines if the point p is inside the capsule C
%
% Syntax:
%    result = containsPoint(C,p)
%
% Inputs:
%    C - capsule
%    p - point specified as a vector
%
% Outputs:
%    res - true/false
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    p = [1; 1; 1];
%    res = contains(C,p)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-March-2019
% Last update:   15-September-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

numPoints = size(p,2);
% set result to 0
res = false(numPoints,1);

for iPoint = 1:numPoints
    p_curr = p(:,iPoint);
    
    % in case capsule is just a ball
    if isempty(C.g)
        p_delta = p_curr-C.c; % point relative to center
        if norm(p_delta) <= C.r
            res(iPoint) = 1; % p is contained in ball
        end
    % capsule is not a ball
    else
        % normalized generator
        g_length = norm(C.g);
        g_norm = C.g/g_length;

        % check if projection of p lies on axis of hyper-cylinder
        p_delta = p_curr-C.c; % point relative to center
        p_dir = p_delta'*g_norm;
        p_proj = p_dir*g_norm;

        % does p lie on axis of hyper-cylinder?
        if abs(p_dir) <= g_length 
            % is perpendicular part smaller than the radius?
            p_perpendicular = p_delta - p_proj;
            if norm(p_perpendicular) <= C.r
                res(iPoint) = true; % p is contained in cylinder and thus in  capsule
            end
        % is p in up part of both half-hyperballs?
        elseif p_dir > 0
            p_delta_up = p_curr-(C.c + C.g);
            if norm(p_delta_up) <= C.r
                res(iPoint) = true; % p is contained in up part of half-hyperballs
            end
        % is p in down part of both half-hyperballs?
        elseif p_dir < 0
            p_delta_dw = p_curr-(C.c - C.g);
            if norm(p_delta_dw) <= C.r
                res(iPoint) = true; % p is contained in down part of half-hyperballs
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
