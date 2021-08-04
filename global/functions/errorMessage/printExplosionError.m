function ME = printExplosionError()
% printExplosionError - standardized MException object if reachable sets
%    explode during the analysis; this function is necessary to ensure
%    that all strcmp() match
%
% Syntax:  
%    printExplosionError()
%
% Inputs:
%    -
%
% Outputs:
%    ME - MException object
%
% Example: 

% Author:       Mark Wetzlinger
% Written:      19-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

ME = MException('reach:setexplosion',...
                'Abort analysis due to reachable set explosion!');

end

