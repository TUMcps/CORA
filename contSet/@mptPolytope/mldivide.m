function obj = mldivide(obj,P)
% mldivide - computes the set difference obj\P such that P is subtracted
% from obj
%
% Syntax:  
%    obj = mldivide(obj,P)
%
% Inputs:
%    obj - mptPolytope object
%    P - mptPolytope object
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
% Written:      03-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute set difference
obj.P = obj.P\P.P;


%------------- END OF CODE --------------