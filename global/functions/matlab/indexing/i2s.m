function subscriptMatrix = i2s(siz,indexVector)
% i2s - extends the MATLAB function ind2sub so that it can be used
% conveniently for arbitrarily many dimensions. The function ind2sub 
% converts linear indices to subscripts
%
% Syntax:
%    subscriptMatrix = i2s(siz,indexVector)
%
% Inputs:
%    siz - vector of number of segments for each dimension as a row vector
%    ("size")
%    indexVector - vector of linear indices that should be converted into
%    subscripts
%
% Outputs:
%    subscriptMatrix - matrix of subscripts; each row corresponds to one
%    index
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       14-September-2006
% Last update:   28-July-2020
%                30-September-2024 (TL, rewrote)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extend given size to have at least two elements
if isscalar(siz)
    dims = [siz,1];
else
    dims = siz;
end

% compute subscripts
subs = cell(numel(dims),1);
[subs{:}] = ind2sub(dims,indexVector);

% rearrange such that each row corresponds to one index pair
subscriptMatrix = cell2mat(subs)';

% remove additional columns added due to scalar input
subscriptMatrix = subscriptMatrix(:,1:numel(siz));
    
end
    
% ------------------------------ END OF CODE ------------------------------
