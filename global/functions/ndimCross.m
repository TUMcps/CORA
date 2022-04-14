function [v] = ndimCross(Q)
% ndimCross - Computes the n-dimensional cross product
%
% Syntax:  
%    [v] = nDimCross(Q)
%
% Inputs:
%    Q - 
%
% Outputs:
%    v - 
%
% Example: 
%    -
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

for i=1:length(Q(:,1))
    D=Q;
    D(i,:)=[];
    v(i,1)=(-1)^(i+1)*det(D);
end

%------------- END OF CODE --------------