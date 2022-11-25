function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if halfspace obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - halfspace object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    hs = halfspace([-1;-1],0);
%    zono1 = zonotope([0 1 1 0; 0 1 0 1]);
%    zono2 = zono1 - [3;3];
%
%    isIntersecting(hs,zono1)
%    isIntersecting(hs,zono2)
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
% See also: conHyperplane/isIntersecing

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  14-Sep-2019
%               20-Nov-2019
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    % get halfspace object
    if ~isa(obj1,'halfspace')
       temp = obj1;
       obj1 = obj2;
       obj2 = temp;
    end
    
    % check user input for correctness
    if strcmp(type,'exact')
       if isa(obj2,'taylm') || isa(obj2,'polyZonotope') || ... 
          isa(obj2,'ellipsoid') || isa(obj2,'capsule') 
      
            error('No exact algorithm implemented for this set representation!');
       end
    end
    
    % check for intersection
    bound = supportFunc(obj2,obj1.c,'lower');
    res = bound <= obj1.d;

%------------- END OF CODE --------------