function [qZsplit] = splitLongestGen(qZ)
% splitLongestGen - Splits the longest generator factor of a quadZonotope
%
% Syntax:  
%    [qZsplit] = splitOneGen(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    qZsplit - cell array of split quadZonotopes
%
% Example: 
%    ---
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%compute metric of generators (shortest generators)
h=vnorm(qZ.G,1,2);
[elements,indices]=sort(h);

%split
qZsplit = splitOneGen(qZ, indices(end));


%------------- END OF CODE --------------