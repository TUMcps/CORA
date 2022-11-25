function [obj] = halfspace(obj)
% halfspace - dummy function for creating halfspace representation
%
% Syntax:  
%    [obj] = halfspace(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    obj - mptPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      12-February-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve halfspace representation
try %MPT V3
    Hfull = obj.P.H;
    H = Hfull(:,1:end-1);
    K = Hfull(:,end);
catch %MPT V2
    [H,K] = double(obj.P);
end
%write to object structure
obj.halfspace.H=H;
obj.halfspace.K=K;
obj.halfspace.equations=length(K);

%------------- END OF CODE --------------