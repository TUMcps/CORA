function obj = and(obj,S)
% and - computes the intersection of a constrained hyperplane with a set
%
% Syntax:  
%    obj = and(obj,S)
%
% Inputs:
%    obj - conHyperplane object
%    S - conSet object
%
% Outputs:
%    obj - conSet object
%
% Example: 
%    poly = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    cH = conHyperplane([1 1],2,[-1 0],-1);
%
%    res = cH & poly;
%
%    figure
%    hold on
%    xlim([-2,4]);
%    ylim([-4,4]);
%    plot(cH,[1,2],'r','LineWidth',3);
%    plot(poly,[1,2],'b');
%    plot(res,[1,2],'g');
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