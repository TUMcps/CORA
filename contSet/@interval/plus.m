function res = plus(summand1,summand2)
% plus - Overloaded '+' operator for intervals
%
% Syntax:  
%    res = plus(summand1,summand2)
%
% Inputs:
%    summand1 - interval (for computational efficiency, no single value
%    considered; does not require type checking)
%    summand2 - interval (for computational efficiency, no single value
%    considered; does not require type checking)
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      19-June-2015
% Last update:  23-June-2015
%               10-August-2016
%               24-August-2016
% Last revision:---

%------------- BEGIN CODE --------------


%Find an interval object
%Is summand1 an interval?
if isa(summand1,'interval')
    
    %initialize resulting interval
    res = summand1;
    %initialize other summand
    summand = summand2;
    
%Is summand2 an interval?    
elseif isa(summand2,'interval')
    
    %initialize resulting zonotope
    res = summand2;
    %initialize other summand
    summand = summand1;
       
end

%%Is summand an interval?
if isa(summand,'interval')
    %Calculate infimum and supremum
    res.inf = res.inf + summand.inf;
    res.sup = res.sup + summand.sup;

% is summand a zonotope
elseif isa(summand,'zonotope')
    res = zonotope(res) + summand;
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate infimum and supremum
    res.inf = res.inf + summand;
    res.sup = res.sup + summand;
    
%something else?    
else
    res.inf=[];
    res.sup=[];
    error('this operation is not implemented');
end

%------------- END OF CODE --------------