function [yjump] = reset(obj,y0)
% reset - resets the continuous state according the linear reset function
%
% Syntax:  
%    [yjump] = reset(obj,y0)
%
% Inputs:
%    obj - transition object
%    y0 - state value before the reset
%
% Outputs:
%    yjump - state value after the reset
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 04-May-2007 
% Last update: 07-October-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from object structure
A=obj.reset.A;
b=obj.reset.b;

%get correct orientation of y0
[n,m]=size(y0);
if m~=1
    y0=y0';
end

%if y0 is not a cell
if ~iscell(y0)
    yjump=A*y0+b;
else
    for i=1:length(y0)
        if ~isempty(y0{i})
            yjump{i}=A*y0{i}+b;
        end
    end
end

%------------- END OF CODE --------------