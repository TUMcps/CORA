function [yjump] = resetTrans(obj,y0)
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

% Author:       Matthias Althoff
% Written:      13-August-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object structure
b=obj.reset.b;

%get correct orientation of y0
[n,m]=size(y0);
if m~=1
    y0=y0';
end

%if y0 is not a cell
if ~iscell(y0)
    yjump=y0+b;
else
    for i=1:length(y0)
        yjump{i}=y0{i}+b;
    end
end

%------------- END OF CODE --------------