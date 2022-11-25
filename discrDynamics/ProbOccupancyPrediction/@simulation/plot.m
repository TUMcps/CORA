function plot(Obj)
% Purpose:  Play movie of space probabilities 
% Pre:      Simulation object
% Post:     ---

%get field
field=get(Obj.markovchain,'field');
%generate rectangles-----------------------------------------
for cellNr=1:(length(Obj.probabilityDistribution(:,1))-1)
    %generate polytope out of cell
    [IHP]=segmentPolytope(field,cellNr); %IHP:interval hull polytope 
    V=extreme(IHP);
    xV=V(:,1);
    yV=V(:,2);
    k=convhull(xV,yV);
    
    x{cellNr}=xV(k);
    y{cellNr}=yV(k);
end
colormap('pink')
for k=1:length(Obj.probabilityDistribution(1,:))
    for i=1:(length(Obj.probabilityDistribution(:,1))-1)
        fill(x{i},y{i},Obj.probabilityDistribution(i+1,k));
        hold on
    end
    M(k)=getframe;
end

movie(M,1,5);
movie(M,1,5);
