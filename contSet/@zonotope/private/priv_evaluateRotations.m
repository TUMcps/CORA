function [iAxis,relImpr,direction] = priv_evaluateRotations(V,deltaAngle)
% priv_evaluateRotations - check which rotation decreases the volume of an
%    enclosing bounding box
%
% Syntax:
%    [iAxis,relImpr,direction] = priv_evaluateRotations(V,deltaAngle)
%
% Inputs:
%    V - matrix of points
%    deltaAngle - increment value of the angles
%    direction - ???
%
% Outputs:
%    iAxis - indicates which axis has been rotated
%    relImpr - computes the relative improvement of the rotation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-October-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute volume before rotation
I = interval.enclosePoints(V);
n = dim(I);
origVol = volume_(I);

% %obtain rotation center
% c=center(IH);
% %translate the center to the origin
% V=V-c;

%perform rotations
vol1 = zeros(n-1,1);
vol2 = zeros(n-1,1);
for iRot=1:(n-1)
    % obtain projected dimensions
    dims = [iRot iRot+1];
    % rotate vertices
    Vrot1 = rotate(V,dims,deltaAngle);
    Vrot2 = rotate(V,dims,-deltaAngle);
    
    % compute volume of enclosing interval
    vol1(iRot) = volume_(interval(Vrot1));
    vol2(iRot) = volume_(interval(Vrot2));
end

[val1,indx1] = sort(vol1);
[val2,indx2] = sort(vol2);

if val1(1) < val2(1)
    direction=1;
    iAxis=indx1(1);
    relImpr=(origVol-vol1(iAxis))/origVol;
else
    direction=-1;
    iAxis=indx2(1);
    relImpr=(origVol-vol2(iAxis))/origVol;
end

% ------------------------------ END OF CODE ------------------------------
