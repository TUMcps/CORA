function [mat]=vec2mat(varargin)
% vec2mat - Stores entries of a vector in a matrix
%
% Syntax:  
%    [mat]=vec2mat(varargin)
%
% Inputs:
%    vec - vector
%    cols - number of columns for the matrix
%
% Outputs:
%    mat - matrix
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
% Last update:  22-June-2010
%               05-October-2010
% Last revision:---

%------------- BEGIN CODE -------------


if nargin==1
    vec = varargin{1};
    %get columns
    cols = sqrt(length(vec));
elseif nargin==2
    vec = varargin{1};
    %get columns
    cols = varargin{2};
end

%seperate columns
rows=length(vec)/cols;
for i=1:cols
    mat(:,i)=vec(1:rows);
    vec(1:cols)=[];
end


%------------- END OF CODE --------------