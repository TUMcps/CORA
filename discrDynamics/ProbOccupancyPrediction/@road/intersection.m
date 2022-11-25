function [p]=intersection(R1,R2,p1,p2,devProb1,devProb2,fArray,type)
% intersection - computes the probability that two reachable sets 
% of traffic participants intersect.
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

% Author: Matthias Althoff
% Written:      13-March-2008 
% Last update:  01-April-2008
%               22-April-2008
%               17-June-2009
%               18-August-2016
%               14-August-2018
% Last revision:---

%------------- BEGIN CODE --------------

% figure;
% hold on

if strcmp(type,'car') || strcmp(type,'carCons')
    %obtain radii of enclosing circles of road segments
    r=fArray.radius.r.car;
    rComb=fArray.radius.rComb.car;
    %load intersection database
    intDatabase=fArray.val.car;
    if strcmp(type,'carCons') % choose between normal and conservative intersection
        intFlag=1;
    else 
        intFlag=0;
    end
elseif strcmp(type,'bicycle')
    %obtain radii of enclosing circles of road segments
    r=0.5*(fArray.radius.r.car+fArray.radius.r.bicycle);
    rComb=0.5*(fArray.radius.rComb.car+fArray.radius.rComb.bicycle);
    %load intersection database
    intDatabase=fArray.val.bicycle;    
    if strcmp(type,'carCons') % choose between normal and conservative intersection
        intFlag=1;
    else 
        intFlag=0;
    end
end    

%compute segment center vectors of R1
ind1=find(p1>0);
[c1]=segCenter(R1,ind1);

%compute segment center vectors of R2
ind2=find(p2>0);
[c2]=segCenter(R2,ind2);


%init potPairs
potPairs=[];


%obtain potential segment intersection pairs
if (~isempty(c1) && ~isempty(c2))
    for i=1:length(c1(:,1))
        for j=1:length(c2(:,1))
            dist=norm(c1(i,:)-c2(j,:));
            if dist<rComb
                potPairs(end+1,:)=[ind1(i),ind2(j)];
            end
        end
    end
end


%init probability
p=0;


if ~isempty(potPairs)
    for iPair=1:length(potPairs(:,1))
        %compute deviation numbers whose polytopes might intersect
        [potDevPairs]=devNumbers(R1,R2,potPairs(iPair,:),r,devProb1,devProb2);
        
        %potPairs(iPair,:)
        %potDevPairs
        
        if ~isempty(potDevPairs)
            for k=1:length(potDevPairs(:,1))
                %generate polytopes
                [pos1,angle1]=segData(R1,potPairs(iPair,1),potDevPairs(k,1));
                [pos2,angle2]=segData(R2,potPairs(iPair,2),potDevPairs(k,2));

                %obtain intersection probability
                intersected = intersection_database(R1,fArray,intDatabase,pos1,angle1,pos2,angle2);     

                if intersected>0
                    %compute partial probabilities
                    prob1=p1(potPairs(iPair,1))*devProb1(potDevPairs(k,1));
                    prob2=p2(potPairs(iPair,2))*devProb2(potDevPairs(k,2));

                    %add to probability
                    if intFlag
                        p=p+prob1*prob2;
                    else
                        p=p+intersected*prob1*prob2;
                    end

                    %plot both rectangles for a test
                    %testIntersectionProbability(R1,R2,potPairs(iPair,1),potPairs(iPair,2),...
                    %    potDevPairs(k,1),potDevPairs(k,2),zonotope(fArray.IH.car),Zvehicle,intersected)
                end
            end
        end
    end
end


%-----------------------------------------------------------------
%auxiliary functions
function [potPairs]=devNumbers(R1,R2,potPair,r,devProb1,devProb2)

%get deviation segments with nonzero probability
ind1=find(devProb1);
ind2=find(devProb2);

%compute centers of deviation polytopes
[cDev1]=segCenter(R1,potPair(1),ind1);
[cDev2]=segCenter(R2,potPair(2),ind2);

%init potPairs
potPairs=[];

%obtain potential deviation polytope intersection pairs
for i=1:length(cDev1(:,1))
    for j=1:length(cDev2(:,1))
        dist=norm(cDev1(i,:)-cDev2(j,:));
        if dist<r
            potPairs(end+1,:)=[ind1(i),ind2(j)];
        end
    end
end


%-----------------------------------------------------------------
%test functions

function testTransformation(Zvehicle,angle1,angle2,relAngle,pos1,pos2,relPos)

%generate zonotopes from original coordinates/angles
Z1orig=rotMatrix(angle1)*Zvehicle+pos1;
Z2orig=rotMatrix(angle2)*Zvehicle+pos2;

Z1trans=Zvehicle;
Z2trans=rotMatrix(relAngle)*Zvehicle+relPos;

%plot original constellation
h1=figure;
subplot(1,2,1)
plot(Z1orig);
plot(Z2orig);
axis equal

%plot transformed constellation
subplot(1,2,2)
plot(Z1trans);
plot(Z2trans);
axis equal

close(h1);


function mat=rotMatrix(angle)
mat=[cos(angle) -sin(angle);...
     sin(angle) cos(angle)];

function testIntersectionProbability(R1,R2,iSeg1,iSeg2,iDev1,iDev2,Z1,Z2,prob)

h=figure;
plotVehSeg(R1,iSeg1,iDev1,Z1);
plotVehSeg(R2,iSeg2,iDev2,Z2);
v=axis; 
text(v(1)+0.2*(v(2)-v(1)),v(3)+0.8*(v(4)-v(3)),['p=',num2str(prob,'%1.3f')])
figure(h);
close(h);


%------------- END OF CODE --------------