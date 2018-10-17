function [Zbundle] = enclosePolytopes(obj,P,options)
% enclosePolytopes - encloses a set of polytopes using differnt
% overapproximating zonotopes.
%
% Syntax:  
%    [Zencl] = enclosePolytopes(obj,R0,options)
%
% Inputs:
%    obj - location object
%    P - initial reachable polytopes
%    options - options struct
%
% Outputs:
%    Zencl - cell array of enclosing zonotopes
%
% Example: 
%
% Other m-files required: not specified
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      11-October-2008
% Last update:  31-July-2009
%               27-July-2016
%               24-November-2017
% Last revision:---

%------------- BEGIN CODE --------------


%obtain new vertices
V=[];
timeSteps=length(P);
for iStep=1:timeSteps
    if ~iscell(P{iStep}) % no split sets
        Vpartial=vertices(P{iStep});
        Vcombined=get(Vpartial,'V');
        V=[V,Vcombined];
    else % sets are split
        Vcombined=[];
        for iSubSet=1:length(P{iStep})
            Vpartial=vertices(P{iStep}{iSubSet});
            Vcombined=[Vcombined,get(Vpartial,'V')];
        end
        V=[V,Vcombined];
    end
    if iStep==1
        x0A=mid(interval(vertices(Vcombined)));
    end
    if iStep==timeSteps
        x0B=mid(interval(vertices(Vcombined)));
    end
end

%get previous direction
dim=length(V(:,1));
if ~isempty(x0A) && isempty(x0B)
    dif=x0B-x0A;
    if norm(dif)==0
        prevDirection=[1;zeros(dim-1,1)];
    else
        prevDirection=dif/norm(dif);
    end
else
    prevDirection=[1;zeros(dim-1,1)];
end
    
V=vertices(V);

nrOfEncl=length(options.enclosureEnables);
for i=1:nrOfEncl
    %perform different enclosures
    switch options.enclosureEnables(i)
        case 1
            %obtain direction of the vector field
            x0=mid(interval(V));
            fcnHandle=getfcn(obj.contDynamics,options);
            direction=fcnHandle(0,x0); 
            %new direction
            Ztmp = dirPolytopeNew(V,direction);
        case 2
            %old direction
            Ztmp = dirPolytopeNew(V,prevDirection); 
        case 3
            %box enclosure
            Ztmp = zonotope(interval(V));
        case 4
            %exponential matrix transformation
            Ztmp = options.W{2}*zonotope(interval(pinv(options.W{2})*V));
            pinv(options.W{2})
        case 5
            %principal component analysis
            %Ztmp = options.W{1}*zonotope(pinv(options.W{1})*V);
            Ztmp = zonotope(V);
        otherwise
            disp('No proper enclosure method selected.')
    end
    %store result
    Zencl{i}=Ztmp;
end


% %split first zonotope
% Znew=split(Zencl{1},1);
% %copy other zonotopes
% Zencl{1}=Znew{1};
% Zencl{nrOfEncl+1}=Znew{2};
% for i=2:nrOfEncl
%     Zencl{nrOfEncl+i}=Zencl{i};
% end
%     
% %convert to Zbundle
% %Zbundle=zonotopeBundle(Zencl(1:4));
% Zbundle{1}=zonotopeBundle(Zencl(1:4));
% Zbundle{2}=zonotopeBundle(Zencl(5:8));

if length(Zencl) == 1
    Zbundle = Zencl{1};
else 
    Zbundle=zonotopeBundle(Zencl(:));
end

%test
if isfield(options,'debug') && options.debug==1
    figure;
    hold on
    dims=[1 2];
    plot(Zbundle,dims,'b');
    for iStep=1:length(P)
        plot(P{iStep},dims,'g');
    end
end

% 
% figure;
% hold on
% dims=[2 3];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% figure;
% hold on
% dims=[4 5];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end

%test
if isfield(options,'debug') && options.debug==1
    figure;
    hold on
    dims=[1 2];
    for i=1:length(Zbundle.Z)
        plot(Zbundle.Z{1},dims,'b');
    end
    for iStep=1:length(P)
        plot(P{iStep},dims,'g');
    end
end
% 
% figure;
% hold on
% dims=[2 3];
% plot(Zbundle{1},dims);
% plot(Zbundle{2},dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% figure;
% hold on
% dims=[4 5];
% plot(Zbundle{1},dims);
% plot(Zbundle{2},dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end


% 
% end

% % if options.t>2
% % 
% figure;
% hold on
% dims=[1 2];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% plot(Zencl{4},dims);
% p1=0*x0;
% p2=0*x0;
% p1(dims(1))=1;
% p2(dims(2))=1;
% plot([p1';p2']*V,'r');
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% 
% figure;
% hold on
% dims=[2 3];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% plot(Zencl{4},dims);
% p1=0*x0;
% p2=0*x0;
% p1(dims(1))=1;
% p2(dims(2))=1;
% plot([p1';p2']*V,'r');
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% 
% figure;
% hold on
% dims=[4 5];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% plot(Zencl{4},dims);
% p1=0*x0;
% p2=0*x0;
% p1(dims(1))=1;
% p2(dims(2))=1;
% plot([p1';p2']*V,'r');
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
% % 
% % end

%------------- END OF CODE --------------