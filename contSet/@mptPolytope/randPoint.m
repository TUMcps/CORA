function p = randPoint(P)
% randPoint - generates a random point inside a polytope
%
% Syntax:  
%    P = randPoint(P)
%
% Inputs:
%    p - mptPolytope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    poly = mptPolytope.generateRandom(2);
%
%    points = zeros(2,100);
%    for i = 1:100
%       points(:,i) = randPoint(poly);
%    end
%
%    figure; hold on;
%    plot(poly,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPointExtreme

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % draw n+1 random extreme points
    n = dim(P);
    points = zeros(n,n+1);
    
    for i = 1:size(points,2)
       points(:,i) = randPointExtreme(P); 
    end

    % interpolate betwenn the points
    delta = rand(n+1,1);
    delta = delta./sum(delta);
    
    p = points*delta;

end

%------------- END OF CODE --------------