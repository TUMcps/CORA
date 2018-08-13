function [result] = in(obj1,obj2)
% in - determines if elements of a zonotope obj2 are in another zonotope
% obj1
%
% Syntax:  
%    [result] = in(obj1,obj2)
%
% Inputs:
%    obj1 - 1st zonotope object
%    obj2 - 2nd zonotope object
%
% Outputs:
%    result - 1/0 if zonotope is in, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%simple test: Is the center of obj2 in obj1?
c=center(obj2);
inequality=(obj1.halfspace.H*c<=obj1.halfspace.K);

if all(inequality)
    result=1;
else
end

%use idea of the previous interval test
%make a polytope of obj1
P=polytope(obj1.halfspace.H,obj1.halfspace.K);

%test 1: Is the center of obj2 in obj1?
c=center(obj2);
bool=isinside(P,c);

if bool==1
    result=1;
else
    %test 2: Is any vertex of obj2 in obj1?
    %overapproximate zonotope to order 1 to prevent an explosion of the
    %number of vertices
    obj2=reduce(obj2,'girard',1);
    %get vertices of the zootope
    V=vertices(obj2);
    %check if vertex in polytope
    bool=0;
    iVertex=1;
    while (bool==0) & (iVertex<=length(V(:,1)))
        bool=isinside(P,V(iVertex,:));
        iVertex=iVertex+1;
    end
    if bool==1
        result=1;
    else
        %test 3: check if the intersection of obj1 and obj2 is empty
        intersection=P&polytope(obj2);
        if ~isempty(intersection)
            result=1;
        else
            result=0;
        end
    end
end


%------------- END OF CODE --------------