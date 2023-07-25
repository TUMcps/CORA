function obj = and_(obj,P,varargin)
% and_ - computes intersection of two mptPolytopes
%
% Syntax:  
%    obj = and_(obj,P)
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
% See also: conZonotope/and_

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  ---
% Last revision:27-March-2023 (MW, rename and_)

%------------- BEGIN CODE --------------

% convert second object to mptPolytope
if isa(P,'levelSet')
    obj = and_(P,obj,'exact');

elseif isa(P,'mptPolytope')
    obj.P = obj.P & P.P;

elseif isa(P,'fullspace')
    % keep object as is...

elseif isa(P,'emptySet')
    % intersection is empty set
    obj = P;

else
    % conversion to polytope, then compute intersection
    P = mptPolytope(P);
    obj.P = obj.P & P.P;

end

%------------- END OF CODE --------------