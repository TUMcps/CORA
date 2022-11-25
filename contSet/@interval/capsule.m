function C = capsule(obj)
% capsule - Encloses an interval with a capsule
%
% Syntax:  
%    C = capsule(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    C - capsule object
%
% Example: 
%    I = interval.generateRandom(2);
%    C = capsule(I);
%
%    figure
%    hold on
%    plot(I,[1,2],'r');
%    plot(C,[1,2],'b');
%    axis equal
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:        Niklas Kochdumper
% Written:       23-December-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
    
% dimension with largest width -> generator of capsule
width = rad(obj);
[~,ind] = sort(width,'descend');

g = zeros(dim(obj),1);
g(ind(1)) = width(ind(1));

% radius of capsule = enclosing radius of remaining interval
int_ = project(obj,ind(2:end));
r = radius(int_);

% construct capsule object
C = capsule(center(obj),g,r);

%------------- END OF CODE --------------