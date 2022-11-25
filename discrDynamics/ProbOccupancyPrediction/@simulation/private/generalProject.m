function [P]=generalProject(nrOfSegments,projected_coordinates)
% Purpose:  build projection matrix
% Pre:      nrOfSegments, projected coordinares
% Post:     projection matrix


% Get total number of segments----------------------
totalNrOfSegments=prod(nrOfSegments);
%---------------------------------------------------

%build segment index matrix-------------------------
for i=1:totalNrOfSegments
    list(i,:)=i2s(nrOfSegments',i);
end
%---------------------------------------------------

%get new number of segments
new_nrOfSegments=nrOfSegments(projected_coordinates);

%initialize P---------------------------------------
P(prod(new_nrOfSegments),totalNrOfSegments)=0;
%---------------------------------------------------

%Generate new subscripts for the remaining variables------------------
list2=i2s(new_nrOfSegments',1:prod(new_nrOfSegments));
%---------------------------------------------------------------------

%get indices of projected cells
for iList=1:length(list2(:,1))
    indices1=find(list2(iList,1)==list(:,1));
    indices2=find(list2(iList,2)==list(indices1,2));
    indices=indices1(indices2);
    P(iList,indices)=1;
end

%add outside area aspect
P=sparse([1,0*P(1,:);0*P(:,1),P]); 