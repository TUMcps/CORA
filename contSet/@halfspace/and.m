function obj = and(obj,S)
% and - computes the intersection of a halfspace with a set
%
% Syntax:  
%    obj = and(obj,S)
%
% Inputs:
%    obj - halfspace object
%    S - conSet object
%
% Outputs:
%    obj - conSet object
%
% Example: 
%    % define sets
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    hs = halfspace([1 1],2);
%
%    % intersection
%    res1 = hs & zB;
%
%    % visualization
%    figure
%    hold on
%    xlim([-1,4]);
%    ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res1,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(zB,[1,2],'b','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/and

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    obj = S & obj;

%------------- END OF CODE --------------