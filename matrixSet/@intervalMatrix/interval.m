function I = interval(intMat)
% interval - Converts an interval matrix to an interval vector
%
% Syntax:  
%    I = interval(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    I - interval object
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      21-June-2010
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%convert matrix limits
leftLimit = mat2vec(infimum(intMat.int));
rightLimit = mat2vec(supremum(intMat.int));
    
%instantiate interval hull
I=interval(leftLimit,rightLimit);

%------------- END OF CODE --------------