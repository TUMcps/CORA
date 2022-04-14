function obj = and(obj,P)
% and - computes intersection of two mptPolytopes
%
% Syntax:  
%    obj = and(obj,P)
%
% Inputs:
%    obj - mptPolytope object
%    P - mptPolytope object
%
% Outputs:
%    obj - mptPolytope object
%
% Example: 
%    poly = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    cH = conHyperplane([1 1],2,[-1 0],-1);
%
%    res = poly & cH;
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

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % convert second object to mptPolytope
    if ~isa(P,'mptPolytope')
        if isa(P,'levelSet')
            obj = P & obj;
            return;
        else
            P = mptPolytope(P); 
        end
    end

    % compute intersection
    obj.P = obj.P & P.P;

%------------- END OF CODE --------------