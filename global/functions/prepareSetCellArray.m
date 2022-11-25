function [S,size_S] = prepareSetCellArray(S,obj)
% prepareSetCellArray - Checks validity of S against obj (set representation)
%
% Syntax:  
%    S=prepareSetCellArray(S,obj)
%
% Inputs:
%    S - array (or single instance) of "stuff"
%    obj - valid set representation
%
% Outputs:
%    S- collapsed and checked cell array
%    size_S - original size of S
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Victor Gassmann
% Written:      15-March-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE -------------
if ~ismethod(obj,'dim')
    error('obj is not a valid set representation');
end
size_S = [0,0];
if ~iscell(S)
    if isempty(S)
        return;
    end
    if isa(S,'double')
        S = mat2cell(S,size(S,1),ones(1,size(S,2)));
    else
        S = {S};
    end
end

if isempty(S)
    return;
end
size_S = size(S);
S = S(:);

n = dim(obj);

% check if all elements of S are of the same type
if length(S)>1 && ~all(cellfun(@(ss)isa(ss,class(S{1})),S))
    error('All elements of cell array need to be of same type!');
end
% check dimensions
valid = true;
for i=1:length(S)
    if isa(S{i},'double')
        if length(S{i})~=n && length(S{i})==1
            S{i} = repmat(S{i},n,1);
        end
        valid = valid && length(S{i}) == n;
    else
        valid = valid && dim(S{i}) == n;
    end
    if ~valid
        error('Dimensions of obj and elements of array do not match!');
    end
end

%------------- END OF CODE --------------