function indices = s2i(siz,subscripts)
% s2i - extends the MATLAB function sub2ind so that it can be used
%    conveniently for arbitrarily many dimensions. The function sub2ind 
%    converts subscripts to linear inidces
%
% Syntax:
%    indices = s2i(siz,subscripts)
%
% Inputs:
%    siz - vector of number of segments for each dimension as a row vector
%          ("size")
%    subscripts - matrix of subscripts that should be converted into linear 
%                 indices; each row is converted into a different index
%
% Outputs:
%    indices - vector of indices
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       14-September-2006
% Last update:   28-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

multiplicator = siz;
for i = 2:length(multiplicator)
    multiplicator(i) = multiplicator(i)*multiplicator(i-1);
end
multiplicator = [1, multiplicator(1:end-1)];

indices = (multiplicator * (subscripts - 1)')+1;
indices=indices.*(prod((subscripts<=siz).*(subscripts>0),2))';

% ------------------------------ END OF CODE ------------------------------
