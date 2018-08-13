function res = norm(obj, varargin)
% norm - computes exactly the maximum norm value of all points in a
% zonotope
%
% Syntax:  
%    res = norm(obj, varargin)
%
% Inputs:
%    obj - zonotope object
%    varargin - list of optional inputs
%
% Outputs:
%   res - resulting maximum norm value
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
% Written:      29-November-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain center and generators
c = obj.Z(:,1);
G = obj.Z(:,2:end);

% compute point currently providing the maximum norm value
maxPoint = c;

% loop through generators
for iGen = 1:length(G(1,:))
    % check whether addition or subtraction results in larger norm value
    maxPointPlus  = maxPoint + G(:,iGen);
    maxPointMinus = maxPoint - G(:,iGen);
    % plus or minus?
    if norm(maxPointPlus,varargin{:}) > norm(maxPointMinus,varargin{:})
        maxPoint = maxPointPlus;
    else
        maxPoint = maxPointMinus;
    end
end

% the norm of the set equals the norm of the maxPoint
res = norm(maxPoint,varargin{:});


%------------- END OF CODE --------------