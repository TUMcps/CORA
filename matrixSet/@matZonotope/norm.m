function res = norm(obj, varargin)
% norm - computes exactly the maximum norm value of all possible matrices
%
% Syntax:  
%    res = norm(obj, varargin)
%
% Inputs:
%    matZ - matrix zonotope
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
% See also: @zonotope/norm

% Author:       Matthias Althoff
% Written:      02-November-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain center and generators
c = obj.center;
G = obj.generator;

% compute point currently providing the maximum norm value
maxPoint = c;

% loop through generators
failed = 0;
for iGen = 1:length(G)
    % check whether addition or subtraction results in larger norm value
    maxPointPlus  = maxPoint + G{iGen};
    maxPointMinus = maxPoint - G{iGen};
    % plus or minus?
    nPlus = norm(maxPointPlus,varargin{:});
    nMinus = norm(maxPointMinus,varargin{:});
    if nPlus==nMinus
        failed = 1;
        break
    else
        if nPlus > nMinus
            maxPoint = maxPointPlus;
        else
            maxPoint = maxPointMinus;
        end
    end
end

% backup strategy for failed attempt
if failed
    %overapproximate by interval matrix
    M = intervalMatrix(obj);
    
    % obtain absolute values of interval matrix
    M_abs = abs(M);

    % compute norm of absolute value
    res = norm(M_abs, varargin{:});
    
    disp('backup strategy for norm computation of matrix zonotope applied');
else
    % the norm of the set equals the norm of the maxPoint
    res = norm(maxPoint,varargin{:});
end


%------------- END OF CODE --------------