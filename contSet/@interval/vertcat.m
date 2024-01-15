function I = vertcat(varargin)
% vertcat - Overloads the opertor for vertical concatenation
%
% Syntax:
%    I = vertcat(varargin)
%
% Inputs:
%    varargin - list of interval objects 
%
% Outputs:
%    I - interval object 
%
% Example: 
%    I1 = interval(-1, 1);
%    I2 = interval(1, 2);
%    I = [I1,I2];
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       26-June-2015 
% Last update:   08-August-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = varargin{1};

%if object is not an interval
if ~isa(I,'interval')
    tmp = I;
    I = interval(tmp,tmp);
end

for i = 2:nargin
    %check if concatented variable is an interval
    if isa(varargin{i},'interval')
        I.inf = [I.inf; varargin{i}.inf];
        I.sup = [I.sup; varargin{i}.sup];
    else
        I.inf = [I.inf; varargin{i}];
        I.sup = [I.sup; varargin{i}];
    end
end


% ------------------------------ END OF CODE ------------------------------
