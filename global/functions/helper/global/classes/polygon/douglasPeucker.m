function res = douglasPeucker(points,varargin)
% douglasPeucker - simplify a 2D line segment using Douglar-Peucker alg.
%
% Syntax:
%    res = douglasPeucker(points)
%    res = douglasPeucker(points,tol)
%
% Inputs:
%    points - points defining the line segment
%    tol - tolerance for the Hausdorf distance between original and
%          simplified line segment 
%
% Outputs:
%    res - points defining the simplified line segment
%
% Example:
%    % line segment
%    x = [0 1 4 6 9 10 11 15];
%    y = [3 1 0.5 1 2 2 1 3];
%    points = [x;y];
%    
%    % simplify line segment
%    points_ = douglasPeucker(points,0.5);
% 
%    % visualization
%    figure; hold on;
%    plot(points(1,:),points(2,:),'r');
%    plot(points_(1,:),points_(2,:),'--b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon, capsule

% Authors:       Niklas Kochdumper
% Written:       07-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    tol = setDefaultValues({0.01},varargin);

    % find point with the largest distance
    dmax = 0;
    index = 0;
    for i = 2:size(points,2)-1
        d = aux_distance(points(:,i),points(:,1),points(:,end));
        if d > dmax
            index = i;
            dmax = d;
        end
    end

    % maximum distance > epsilon -> simplify with recursive function call
    if dmax > tol
        
        % recursive function call
        recRes1 = douglasPeucker(points(:,1:index), tol);
        recRes2 = douglasPeucker(points(:,index:end), tol);

        % store the results
        res = [recRes1(:,1:end-1),recRes2(:,1:end)];
        
    else
        res = [points(:,1),points(:,end)];
    end
end


% Auxiliary functions -----------------------------------------------------

function d = aux_distance(x,p1,p2)
% distance of the point x from the line segment spanned by p1 and p2
    B = gramSchmidt(p1-p2);
    x_ = B'*(x-p1);
    d = abs(x_(2));
end

% ------------------------------ END OF CODE ------------------------------
