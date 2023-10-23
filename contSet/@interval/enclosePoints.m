function I = enclosePoints(points)
% enclosePoints - enclose a point cloud with an interval
%
% Syntax:
%    I = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    I - interval object
%
% Example: 
%    p = -1 + 2*rand(2,10);
%    I = interval.enclosePoints(p);
%    
%    figure; hold on;
%    plot(I);
%    plot(p(1,:),p(2,:),'.k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Niklas Kochdumper
% Written:       05-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval(min(points,[],2),max(points,[],2));

% ------------------------------ END OF CODE ------------------------------
