function obj = vertcat(varargin)
% vertcat - Overloads the opertor for vertical concatenation, e.g. 
% a = [b;c;d];
%
% Syntax:  
%    obj = horzcat(varargin)
%
% Inputs:
%    varargin - list of interval objects 
%
% Outputs:
%    obj - interval object 
%
% Example: 
%    a=interval(-1, 1);
%    b=interval(1, 2);
%    c = [a,b];
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-June-2015 
% Last update:  08-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

obj = varargin{1};

%if object is not an interval
if ~isa(obj,'interval')
    tmp = obj;
    obj = interval();
    obj.inf = tmp;
    obj.sup = tmp;
end

for i = 2:nargin
    %check if concatented variable is an interval
    if isa(varargin{i},'interval')
        obj.inf = [obj.inf; varargin{i}.inf];
        obj.sup = [obj.sup; varargin{i}.sup];
    else
        obj.inf = [obj.inf; varargin{i}];
        obj.sup = [obj.sup; varargin{i}];
    end
end


%------------- END OF CODE --------------