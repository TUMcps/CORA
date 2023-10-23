function A = emptyClassArray(classname,varargin)
% emptyClassArray - constructs class array of specified dimensions
%
% Syntax:
%    A = emptyClassArray(classname,n)
%    A = emptyClassArray(classname,n1,n2,...)
%    A = emptyClassArray(classname,[n1,n2,...])
%
% Inputs:
%    classname  - class name
%    n          - number of rows and columns
%    ni         - length of i-th dimension
%
% Outputs:
%    A - class array of specified size
%
% Example:
%     A = emptyClassArray('zonotope',5,1,2);

% Authors:       Victor Gassmann
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input arg checks
if isempty(varargin)
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif length(varargin)==1 
    if size(varargin{1},1)~=1
        throw(CORAerror('CORA:wrongValue','first','1 by x size vector'));
    end
    sz = varargin{1};
elseif length(varargin)==2 && any(cellfun(@(cc)numel(cc)~=1,varargin))
    throw(CORAerror('CORA:wrongValue','first/second/third/...','integer'));
else
    sz = horzcat(varargin{:});
end
inputArgsCheck({{classname,'att',{'char','string'},{'istextscalar'}};
                    {sz,'att',{'numeric'},{'integer','nonnegative'}}});

% if one value is empty, use className.empty function to construct
if any(sz==0)
    A = eval(classname).empty(sz);
else
    % convert size to cell array
    sz_c = num2cell(sz);
    A(sz_c{:}) = eval(classname);
end

% ------------------------------ END OF CODE ------------------------------
