function Int = interval(matI)
% interval - Converts an interval matrix to an interval vector
%
% Syntax:  
%    Int = interval(matI)
%
% Inputs:
%    matI - interval matrix
%
% Outputs:
%    Int - interval hull
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-June-2010
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%convert matrix limits
leftLimit = mat2vec(infimum(matI.int));
rightLimit = mat2vec(supremum(matI.int));
    
%instantiate interval hull (using MPT toolbox)
Int=interval(leftLimit,rightLimit);




%------------- END OF CODE --------------