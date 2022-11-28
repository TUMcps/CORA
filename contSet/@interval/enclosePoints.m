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
%    points = -1 + 2*rand(2,10);
%    I = interval.enclosePoints(points);
%    
%    figure; hold on;
%    plot(points(1,:),points(2,:),'.k');
%    plot(int,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:        Niklas Kochdumper
% Written:       05-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

I = interval(min(points,[],2),max(points,[],2));

%------------- END OF CODE --------------