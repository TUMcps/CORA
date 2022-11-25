function set = innerApprox(Z,conv)
% innerApprox - returns an inner-approximation as an object of a given class
%
% Syntax:  
%    set = innerApprox(Z,conv)
%
% Inputs:
%    Z - zonotope object
%    conv - class after conversion
%
% Outputs:
%    set - set after conversion (inner-approximation)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        Mark Wetzlinger
% Written:       10-March-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% dimension and generators
G = generators(Z);
[n,gamma] = size(G);

if strcmp(conv,'interval')
    if n > 2
        error("Interval inner-approximation only for n = 2");
        % method below does not extend to higher dimensions
%         strechFactor = max(sum(abs(G),2));
%         set = enlarge(interval(-ones(n,1),ones(n,1)),strechFactor);
%         return;
    end
    
    % special case for 2D zonotopes
    
    % make first coordinate always positive
    % if first coordinate zero, make second negative
    flipsign = sign(G(1,:)) < 0 | (G(1,:) == 0 & G(2,:) > 0);
    factors = ones(2,gamma) - 2*ones(2,gamma,1) .* flipsign;
    G = G .* factors;

    % order the generators according to y/x ratio
    angle = G(2,:) ./ G(1,:);
    angle(angle == Inf) = -Inf;
    [~,sortIdx] = sort(-angle);
    G = G(:,sortIdx);
    
    % init default bound
    bound = zeros(2,1);
    bestsize = 0;
    
%     figure; hold on; box on;
    
    % loop over all possible paths (always clockwise)
    Gpath = G;
    for tries=1:gamma
        % adapt ordering: first generator times -1, append to end
        if tries > 1
            Gpath = [Gpath(:,2:end), -Gpath(:,1)];
        end
        
        % sum gives end point of path (thus also start point)
        endpoint = sum(Gpath,2);
        startpoint = -endpoint; % due to symmetry
        
        if any(startpoint == 0)
            % if any coordinate is 0, resulting interval would be flat
            continue;
        end
        
        % quadrant of start point
        if all(startpoint > 0)
            quadr = 1; flip = [1;-1];
        elseif all(startpoint < 0)
            quadr = 3; flip = [1;-1];
        elseif startpoint(1) < 0 && startpoint(2) > 0
            quadr = 2; flip = [-1;1];
        else
            quadr = 4; flip = [-1;1];
        end
        % crucial point for containment of box in zonotope
        crucialpoint = flip .* startpoint;
        radcrucialpoint = vecnorm(crucialpoint);
        
        % follow path until know whether crucial point is contained or not
        allvert = startpoint + [zeros(2,1), 2*cumsum(Gpath,2)];
        
        % index where crucialpoint is passed by
        if quadr == 1
            passIdx = find(allvert(1,:) >= crucialpoint(1),1,'last');
        elseif quadr == 2
            passIdx = find(allvert(2,:) >= crucialpoint(2),1,'last');
        elseif quadr == 3
            passIdx = find(allvert(1,:) <= crucialpoint(1),1,'last');
        else
            passIdx = find(allvert(2,:) <= crucialpoint(2),1,'last');
        end
        
        % check containment of crucialpoint in Z
        crucialpointInZ = ~any(passIdx == [1,gamma+1]) ...
            && vecnorm(allvert(:,passIdx)) >= radcrucialpoint ...
            && vecnorm(allvert(:,passIdx+1)) >= radcrucialpoint;
        
%         hold off; plot(Z); hold on;
%         scatter(allvert(1,:),allvert(2,:),16,'r');
%         scatter(crucialpoint(1,:),crucialpoint(2,:),16,'k');
        
        if crucialpointInZ
            thisbound = abs(crucialpoint);
            % take new bound if larger inner-approximation (by area)
            thissize = prod(thisbound);
            if thissize > bestsize
                bound = thisbound;
                bestsize = thissize;
            end
        end
        
    end
    
    % shift by center of zonotope
    bound = bound + center(Z);
    
    % inner-approximation
    set = interval(-bound,bound);
    
else
    % throw error for all other conversions
    error(noops(conv));
end


end

%------------- END OF CODE --------------