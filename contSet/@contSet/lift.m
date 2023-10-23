function S = lift(S,N,proj)
% lift - lifts a set to a higher-dimensional space,
%    having the new dimensions unbounded
%
% Syntax:
%    S = lift(S,N,proj)
%
% Inputs:
%    S - contSet object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional set
%
% Outputs:
%    S - contSet object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/project, contSet/projectHighDim

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 2
    throw(CORAerror("CORA:notEnoughInputArgs",2))
elseif nargin < 3 || isempty(proj)
    proj = 1:dim(S);
end
inputArgsCheck({{S,'att','contSet'};
                {N,'att','numeric',{'nonnan','scalar','nonnegative','integer'}};
                {proj,'att','numeric',{'nonnan','vector','nonnegative'}}});
if dim(S) > N
    throw(CORAerror('CORA:wrongValue','second','Dimension of higher-dimensional space must be larger than or equal to the dimension of the given set.'))
elseif dim(S) ~= length(proj)
    throw(CORAerror('CORA:wrongValue','third','Number of dimensions in higher-dimensional space must match the dimension of the given set.'))
elseif max(proj) > N
    throw(CORAerror('CORA:wrongValue','thrid','Specified dimensions exceed dimension of high-dimensional space.'))
end

% call subfunction
S = lift_(S,N,proj);

% ------------------------------ END OF CODE ------------------------------
