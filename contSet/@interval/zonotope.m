function Z = zonotope(obj)
% zonotope - Converts an interval object into a zonotope object
%
% Syntax:  
%    Z = zonotope(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    I = interval([1;-1], [2; 1]);
%    Z = zonotope(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Matthias Althoff
% Written:      22-July-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % obtain center
    c = center(obj);

    % construct generator matrix G
    G = diag(rad(obj));

    % instantiate zonotope
    Z = zonotope([c,G]);

%------------- END OF CODE --------------