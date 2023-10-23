function [iAxis,relImpr] = evaluateRotations(V,deltaAngle)
% evaluateRotations - check which rotation decreases the volume of an
% enclosing bounding box
%
% Syntax:
%    [iAxis,relImpr] = evaluateRotations(V,deltaAngle)
%
% Inputs:
%    V - matrix of points
%    deltaAngle - increment value of the angles
%
% Outputs:
%    iAxis - indicates which axis has been rotated
%    relImpr - computes the relative improvement of the rotation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       06-October-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute volume before rotation
IH=interval.enclosePoints(V);
origVol=volume_(IH);

% %obtain rotation center
% c=center(IH);
% 
% %translate the center to the origin
% V=V-c;

%get number of dimensions
vertices=V;
dim=length(vertices(:,1));

%perform rotations
for iRot=1:(dim-1)
    %obtain projected dimensions
    dims=[iRot iRot+1];
    %rotate vertices
    Vrot1=rotate(V,dims,deltaAngle);
    Vrot2=rotate(V,dims,-deltaAngle);
    
    %compute enclosing interval
    IH1=interval(Vrot1);
    IH2=interval(Vrot2);
    
    %compute enclosing volume
    vol1(iRot)=volume_(IH1);
    vol2(iRot)=volume_(IH2);
end

[val1,indx1]=sort(vol1);
[val2,indx2]=sort(vol2);

if val1(1)<val2(1)
    direction=1;
    iAxis=indx1(1);
    relImpr=(origVol-vol1(iAxis))/origVol;
else
    direction=-1;
    iAxis=indx2(1);
    relImpr=(origVol-vol2(iAxis))/origVol;
end

% ------------------------------ END OF CODE ------------------------------
