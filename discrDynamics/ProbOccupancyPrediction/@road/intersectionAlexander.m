function [p]=intersectionAlexander(R1,R2,p1,p2,devProb1,devProb2,fArray,deltaV)
% intersection - computes the probability that two reachable sets 
% of traffic participants intersect; for simplicity it is firstly assumed
% that the deviation probability distribution is the same for both traffic
% participants
%
% Syntax:  
%    [p]=intersection(R1,R2,p1,p2)
%
% Inputs:
%    R1 - road object of traffic participant 1
%    R2 - road object of traffic participant 2
%    p1 - segment probability vector of TP 1
%    p2 - segment probability vector of TP 2
%    devProb1 - deviation probability distribution of TP 1
%    devProb2 - deviation probability distribution of TP 2
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

% Author:       Matthias Althoff
% Written:      03-November-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%distinguish between p and pPos!
p1Pos=p1.pos;
p1Total=p1.total;
p2Pos=p2.pos;
p2Total=p2.total;

%load intersection database
intTable=fArray.table;
disc=fArray.disc;

%compute segment center vectors of R1
ind1=find(p1Pos>0);
[c1]=segCenter(R1,ind1);

%compute segment center vectors of R2
ind2=find(p2Pos>0);
[c2]=segCenter(R2,ind2);


%init crash probability
p=0;

%nr of ...
nrOfSegments=length(p1Pos);
nrOfAngleSegments=disc.nrOfAngleSeg;
nrOfxSegments=disc.nrOfxSeg;
nrOfySegments=disc.nrOfySeg;

%velocity segments of the "original system"
nrOfVelocities= (length(p1Total)-1)/nrOfSegments;


for segCounter1=1:length(c1)
    for segCounter2=1:length(c2)
        
        %determine path indices
        iSeg1=ind1(segCounter1);
        iSeg2=ind2(segCounter2);
        
        for iDev1=1:R1.nrOfDevSegments
            for iDev2=1:R2.nrOfDevSegments
       
                %generate polytopes
                [pos1,angle1]=segData(R1,iSeg1,iDev1);
                [pos2,angle2]=segData(R2,iSeg2,iDev2);

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

                %obtain angle, position and velocity segments
                iAngleSeg=round(relAngle/abs(disc.deltaAngle));
                iXseg=round(relPos(1)/abs(disc.deltaX));
                iYseg=round(relPos(2)/abs(disc.deltaY));

                if (iAngleSeg<=nrOfAngleSegments) && (iXseg<=nrOfxSegments) && (iYseg<=nrOfySegments)

                    for iVel1=1:nrOfVelocities
                        for iVel2=1:nrOfVelocities

                            %determine velocity segment of Alexanders table
                            iV1seg=round(deltaV/disc.deltaVA)*iVel1;
                            iV2seg=round(deltaV/disc.deltaVB)*iVel2;

                            %access database as:
                            %Table{iAngleSeg}(iAVseg+NAV*iBVseg+NAV*NBV*inrOfVelocitiesYseg+NAV*NBV*nrOfySeg*iXseg)
                            %count segments indices from 0; add 1 in the end;
                            %use round() to determine indices.
                            tableIndex=1+iV1seg+...
                                       nrOfVelocities*iV2seg+...
                                       nrOfVelocities^2*iYseg+...
                                       nrOfVelocities^2*nrOfySegments*iXseg;
                                   

                            %look up if intersected
                            intersected=intTable{iAngleSeg+1}(tableIndex);


                            if intersected>0
                                %compute indices of total probabilities
                                indTotal1=iSeg1+nrOfSegments*(iVel1-1);
                                indTotal2=iSeg2+nrOfSegments*(iVel2-1);

                                %compute partial probabilities
                                prob1=p1Total(indTotal1)*devProb1(iDev1);
                                prob2=p2Total(indTotal2)*devProb2(iDev2);

                                %add to probability
                                p=p+intersected*prob1*prob2;
                            end
                        end
                    end
                end
            end
        end
    end
end




%------------- END OF CODE --------------