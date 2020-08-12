function res = in(obj1,obj2)
% in - determines if obj2 is fully inside a halfspace obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%
% Inputs:
%    obj1 - halfspace object
%    obj2 - contSet object
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    hs = halfspace([-1;-1],0);
%    zono1 = zonotope([3 1 1 0;3 1 0 1]);
%    zono2 = zono1 - [3;3];
%
%    in(hs,zono1)
%    in(hs,zono2)
%
%    figure
%    hold on
%    xlim([-6,6]);
%    ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(zono1,[1,2],'g','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    xlim([-6,6]);
%    ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(zono2,[1,2],'r','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/in

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-June-2016
% Last update:  27-July-2016
%               02-Sep-2019
%               19-Nov-2019 (NK, extend to all set representations)
% Last revision:---

%------------- BEGIN CODE --------------

    % point in halfspace containment
    if isnumeric(obj2)
        res = obj1.c' * obj2 <= obj1.d;
        
    % set in halfspace containment
    else
        val = supportFunc(obj2,obj1.c,'upper');
        res = val <= obj1.d;
    end

%------------- END OF CODE --------------