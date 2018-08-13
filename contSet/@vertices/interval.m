function I = interval(V)
% interval - Encloses a point set by an interval hull
%
% Syntax:  
%    I = interval(V)
%
% Inputs:
%    V - vertices object
%
% Outputs:
%    I - interval object
%
% Example: 
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      27-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%if vertices object is nonempty
if ~isempty(V.V)
    %minimum values
    leftLimit=min(V.V,[],2);

    %maximum values
    rightLimit=max(V.V,[],2);

    %instantiate interval hull
    I = interval(leftLimit,rightLimit);
else
    I = [];
end


%------------- END OF CODE --------------