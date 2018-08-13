function [IHtotal] = intervalhullAccurate(qZ,splits)
% intervalhullAccurate - Overapproximates a quadZonotope by a tight 
% interval hull
%
% Syntax:  
%    [IHtotal] = intervalhullAccurate(qZ,splits)
%
% Inputs:
%    qZ - quadZonotope object
%    splits - number of splits to improve accuracy
%
% Outputs:
%    IHtotal - intervalhull object
%
% Example: 
%    ---
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      20-September-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------


%splitting 
qZsplit{1} = qZ;

for i=1:splits
    qZnew = [];
    for j=1:length(qZsplit)
        res = splitLongestGen(qZsplit{j});
        qZnew{end+1} = res{1};
        qZnew{end+1} = res{2};
    end
    qZsplit = qZnew;
end

%unify interval hulls
for i=1:length(qZsplit)
    IH=interval(zonotope(qZsplit{i}));
    if i==1
        IHtotal = IH;
    else
        IHtotal = IHtotal | IH;
    end
end


%------------- END OF CODE --------------