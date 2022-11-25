function res = or(Int1,Int2)
% or - compute the union of two intervals
%
% Syntax:  
%    res = or(Int1,Int2)
%
% Inputs:
%    Int1 - interval object
%    Int2 - interval object
%
% Outputs:
%    res - union of the two intervals
%
% Example: 
%    int1 = interval([-2;-2],[-1;-1]);
%    int2 = interval([0;0],[2;2]);
%    
%    res = int1 | int2;
%
%    figure
%    hold on
%    plot(int1,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(int2,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(res,[1,2],'r');
%    xlim([-3,3]);
%    ylim([-3,3]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/split

% Author:       Niklas Kochdumper
% Written:      25-July-2019
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% determine the interval object
if ~isa(Int1,'interval')
    temp = Int1;
    Int1 = Int2;
    Int2 = temp;
end

% different cases depending on the class of the summand
if isa(Int2,'interval')
    
    if dim(Int1) ~= dim(Int2)
        [id,msg] = errDimMismatch();
        error(id,msg);
    end

    res = interval(min([Int1.inf,Int2.inf],[],2), ...
                   max([Int1.sup,Int2.sup],[],2));

elseif isnumeric(Int2)

    res = interval(min([Int1.inf,Int2],[],2), ...
                   max([Int1.sup,Int2],[],2));

elseif isa(Int2,'zonotope') || isa(Int2,'conZonotope') || ...
       isa(Int2,'zonoBundle') || isa(Int2,'polyZonotope') || ...
       isa(Int2,'mptPolytope') || isa(Int2,'conPolyZono')

    res = Int2 | Int1;

else
    % throw error for given arguments
    error(noops(Int1,Int2));
end

%------------- END OF CODE --------------