function int = enclosePoints(points)
% enclosePoints - enclose a point cloud with an interval
%
% Syntax:  
%    int = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    int - interval object
%
% Example: 
%    points = -1 + 2*rand(2,10);
%
%    int = interval.enclosePoints(points);
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(int,[1,2],'r');
%
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

    infi = min(points,[],2);
    sup = max(points,[],2);
    
    int = interval(infi,sup);

%------------- END OF CODE --------------