function intMat = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of matrix or an 
% interval matrix with an interval matrix
%
% Syntax:  
%    intMat = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or interval matrix
%    factor2 - numerical matrix or interval matrix
%
% Outputs:
%    intMatZ - interval matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%set multiplication scheme
if isnumeric(factor1)
    if ~isempty(factor2.setting)
        intvalinit(factor2.setting);
    end
else
    if ~isempty(factor1.setting)
        intvalinit(factor1.setting);
    end
end

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix zonotope
    intMat=factor2;
    %compute new intervals
    intMat.int=matrix*intMat.int;
    
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix zonotope
    intMat=factor1;
    %compute new intervals
    intMat.int=intMat.int*matrix;
    
%both factors are interval matrices
else
    %initialize matrix zonotope
    intMat=factor1;
    %compute interval matrix
    intMat.int=factor1.int*factor2.int;
end

%------------- END OF CODE --------------