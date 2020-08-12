function options = set_volApproxMethod(options)
% set_volApproxMethod - set default value for internal settings 
%                       option.volApproxMethod
%
% Syntax:
%    options = set_volApproxMethod(options)
%
% Inputs:
%    options - options for object
%
% Outputs:
%    options - updated options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      02-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isa(options.R0,'polyZonotope')
    if strcmp(options.alg,'poly')
        if isfield(options,'polyZono')
            if isFullDim(options.R0)
                options.polyZono.volApproxMethod = 'interval';
            else
                options.polyZono.volApproxMethod = 'pca';
            end
        end
    end
end
    
    
end

%------------- END OF CODE --------------