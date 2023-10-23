function matV = vertices(intMat)
% vertices - computes the vertices of an interval matrix
%
% Syntax:
%    matV = vertices(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    matV - cell array of matrix vertices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%conversion to an interval
V = vertices(interval(intMat));

%compute vertices
V = unique(V', 'rows')'; %eliminate vectors that occur multiple times

%convert vertices to matrix vertices
matV=cell(length(V(1,:)),1);
for i=1:length(V(1,:))
    matV{i}=vec2mat(V(:,i));
end

% ------------------------------ END OF CODE ------------------------------
