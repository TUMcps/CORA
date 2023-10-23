function p = polygon(C)
% polygon - Under-approximates a two-dimensional capsule by a polygon and 
%    returns its vertices; mainly used for plotting reasons.
%
% Syntax:
%    p = polygon(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    p - ordered set of points of a polygon
%
% Example: 
%    C = capsule([1; 1], [0; 1], 0.5);
%    p = polygon(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-March-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fixed number of points for half circles
nrOfPoints = 500;

% center of first half circle
if ~isempty(C.g)
    c1 = C.c + C.g;
else
    c1 = C.c;
end

% angles
if ~isempty(C.g) && norm(C.g) > 0
    generatorAngle = atan2(C.g(2),C.g(1));
else
    generatorAngle = 0;
end
startAngle1 = generatorAngle - pi/2;
finalAngle1 = generatorAngle + pi/2;

% generate points of first half circle
angle = linspace(startAngle1,finalAngle1,nrOfPoints);
radius = ones(1,nrOfPoints)*C.r;
[x1,y1] = pol2cart(angle,radius);
x1 = x1 + c1(1); % shift by c1
y1 = y1 + c1(2); % shift by c1

% center of second half circle
if ~isempty(C.g)
    c2 = C.c - C.g;
else
    c2 = C.c;
end

% angles (flip previous values)
startAngle2 = finalAngle1;
finalAngle2 = finalAngle1 + pi;

% generate points of first half circle
angle = linspace(startAngle2,finalAngle2,nrOfPoints);
[x2,y2] = pol2cart(angle,radius);
x2 = x2 + c2(1); % shift by c1
y2 = y2 + c2(2); % shift by c1

% concatenate results (start point added to close polygon)
p = [[x1,x2,x1(1)];[y1,y2,y1(1)]];

% ------------------------------ END OF CODE ------------------------------
