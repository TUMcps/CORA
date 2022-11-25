function p = intersection_database(obj,fArray,intDatabase,pos1,angle1,pos2,angle2)
% intersection_database - checks whether two rectangles intersect based on
% a database; an intersection probability is provided for a given
% uncertainty of the center of the rectangles; The database lookup results
% in a small quantization error
%
% Syntax:  
%    [p]=intersection(R1,R2,p1,p2)
%
% Inputs:
%    road - road object
%    intDatabase - intersection database
%    pos1 - position of first object
%    angle1 - orientation of first object
%    pos2 - position of second object
%    angle2 - orientation of second object
%
% Outputs:
%    p - intersection probability
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written:      18-August-2016   
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%nr of segments
nrOfAngleSegments=length(intDatabase);
[nrOfxSegments,nrOfySegments]=size(intDatabase{1});

%generate rotation matrix
rot=[cos(-angle1) -sin(-angle1);...
     sin(-angle1) cos(-angle1)];
%return relAngle
relAngle=angle2-angle1;
relAngle=rem(relAngle,pi);
%return relPosition
relPos=rot*(pos2-pos1);
%generate positive position differences
if prod(relPos)<0
    relAngle=-relAngle;
end
relPos=abs(relPos);
%generate positive angle differences    
if relAngle<0
    relAngle=relAngle+pi;
end   

%obtain angle and position segments
iAngleSeg=round(relAngle/fArray.segLength.angle)+1;
iXseg=round(relPos(1)/fArray.segLength.x)+1;
iYseg=round(relPos(2)/fArray.segLength.y)+1;

%in case some segments are zero
if (iAngleSeg==0) || (iAngleSeg>nrOfAngleSegments)
    p = 0;
    return;
end
if (iXseg==0) || (iXseg>nrOfxSegments)
    p = 0;
    return;
end   
if (iYseg==0) || (iYseg>nrOfySegments) 
    p = 0;
    return;
end  

%look up intersection probability
p = intDatabase{iAngleSeg}(iXseg,iYseg);


  

%------------- END OF CODE --------------