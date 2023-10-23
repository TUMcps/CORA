function C = capsule(I)
% capsule - Encloses an interval with a capsule
%
% Syntax:
%    C = capsule(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    C - capsule object
%
% Example: 
%    I = interval([-2;-1],[3;4]);
%    C = capsule(I);
%
%    figure; hold on;
%    plot(I,[1,2],'r');
%    plot(C,[1,2],'b');
%    axis equal
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Authors:       Niklas Kochdumper
% Written:       23-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% dimension with largest width -> generator of capsule
width = rad(I);
[~,ind] = sort(width,'descend');

g = zeros(dim(I),1);
g(ind(1)) = width(ind(1));

% radius of capsule = enclosing radius of remaining interval
int_ = project(I,ind(2:end));
r = radius(int_);

% construct capsule object
C = capsule(center(I),g,r);

% ------------------------------ END OF CODE ------------------------------
