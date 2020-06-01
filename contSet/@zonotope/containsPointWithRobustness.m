function [robustness, direction] = containsPointWithRobustness(zonotope, p, earlyStop)
% containsPointWithRobustness - Computes the distance of the point to the zonotope boundary.
%
% If the robustness is positive, the point is contained in the zonotope.
% If the robustness is negative, the point is not contained and the robustness gives the distance of the point to the zonotope.
%
% The function uses optimization to approximate the best direction. Hence,
% the result might contain a numerical error.
%
% Syntax:  
%    result = containsPoint(Z1,p)
%
% Inputs:
%    Z1 - zonotope object
%    p - point specified as a vector
%    earlyStop - stops the direction optimization when the first negativ
%                robustness value is found
%
% Outputs:
%    robustness - distance to zonotope boundary, negative if not included
%    direction - the direction of the closed boundary point of the zonotope to p
%
% Example: 
%    z = zonotope([0 1 0;0 0 1])
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

if nargin < 3
    earlyStop = 0;
end

bound = cos(pi/1000/1000); % controls the precision of the result

center = zonotope.Z(:,1);
generators = zonotope.Z(:,2:end);
dim = size(center, 1);

p = p - center; % translate the problem so that the zonotope has the origin as the center

n = p/norm(p,2); % start with the direction of the point
beta = 0.5;
gamma = 0.5;
s = 1;

% Use gradient decent optimization to find the optimal direction
while 1
    signedGenerators = repmat(sign(n'*generators), dim, 1) .* generators;
    generatorSum = sum(signedGenerators, 2);
    df = generatorSum - p;
    f = n'*df;
    % as we want to have a normed vector, it is better to optimize
    % orthogonaly
    orth = df - n'*df*n;
    isAlreadyOptimal = all(orth == 0);
    if isAlreadyOptimal
        n2 = n;
        f2 = f;
        break;
    end
    orth = orth / norm(orth, 2);
    s = s / beta;

    % Use armijo rule
    % as it is expected that the step size will decrease,
    % reuse the step size over iterations
    while 1
        s = beta * s;
        n2 = n-s*orth;
        n2 = n2/norm(n2, 2);
        f2 = sum(abs(n2'*generators)) - n2'*p;
   
        if n2'*n > bound || (earlyStop && f2 < 0)
            break
        end
        if f2 <= f - gamma*s*df'*orth
            break
        end
    end
    if n2'*n > bound || (earlyStop && f2 < 0)
        break
    end
    n = n2;
end

robustness = f2;
direction = n2;

