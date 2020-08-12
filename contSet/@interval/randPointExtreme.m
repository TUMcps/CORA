function [p] = randPointExtreme(obj)
% randPointExtreme - generates a random extreme point of an interval
%
% Syntax:  
%    [p] = randPointExtreme(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    int = interval.generateRandom(2,[],3);
%    
%    p = randPointExtreme(int);
%
%    figure
%    hold on
%    plot(int,[1,2],'r');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Niklas Kochdumper
% Written:      16-April-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

p = randPointExtreme(zonotope(obj));

%------------- END OF CODE --------------