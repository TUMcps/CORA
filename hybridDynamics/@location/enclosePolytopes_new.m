function [Zbundle] = enclosePolytopes_new(obj,P,options)
% enclosePolytopes - encloses a set of polytopes using differnt
% overapproximating zonotopes.
%
% Syntax:  
%    [Zbundle] = enclosePolytopes(obj,R0,options)
%
% Inputs:
%    obj - location object
%    P - initial reachable polytopes
%    options - options struct
%
% Outputs:
%    Zbundle - zonotope bundle
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
%               26-October-2010
%               17-November-2010
% Last revision:---

%------------- BEGIN CODE --------------


% %obtain prev direction
% x0A=center(parallelotope(P{1}));
% x0B=center(parallelotope(P{end}));

%obtain direction of the vector field
x0=center(parallelotope(P{end}));
fcnHandle=getfcn(obj.contDynamics,options);
direction=fcnHandle(0,x0); 


%specify linear maps
Woption{1} = eye(P{1}.dim);
Woption{2} = mapFromDirection(direction);
Woption{3} = options.W{1};
Woption{4} = options.W{2};

%choose selected maps
W = Woption(options.enclosureEnables);

%ensure that all mapping matrices are normalized
for iMat=1:length(W)
    for iCol=1:P{1}.dim
        W{iMat}(:,iCol)=W{iMat}(:,iCol)/norm(W{iMat}(:,iCol));
    end
end


%initialize enclosing zonotopes
for j=1:length(W)
    Zencl{j}=[];
end

%enclose polytopes
timeSteps=length(P);
for j=1:length(W)
    %initialize centers, genFactors
    centers=[];
    normMax=zeros(1,P{1}.dim);
    for iStep=1:timeSteps
        %enclose single polytopes by zonotopes
        Zpartial=parallelotope(P{iStep},W{j});
        
        %get centers and generators
        Zmat = get(Zpartial,'Z');
        centers(:,end+1)=Zmat(:,1);
        normMax=max(normMax,vnorm(Zmat(:,2:end),1,2));
    end
    %enclose zonotope centers by svd method
    if length(centers(1,:))>1
        Zcenter = zonotope(vertices(centers));
        ZcenterMat = get(Zcenter,'Z');
    else
        ZcenterMat = centers;
    end
    %obtain enclosing generator
    for iCol=1:P{1}.dim
        Gmax(:,iCol)=normMax(iCol)*W{j}(:,iCol);
    end
    %generate enclosing zonotope
    Zencl{j} = zonotope([ZcenterMat,Gmax]);
end


%save result as zonotope bundle
Zbundle=zonotopeBundle(Zencl(:));


% figure;
% hold on
% dims=[1 2];
% plot(Zbundle,dims);
% for iStep=1:length(P)
%     plot(P{iStep},dims,'g');
% end
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

% figure;
% hold on
% dims=[1 2];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% %plot(Zencl{4},dims);
% for iStep=1:timeSteps
%     plot(P{iStep},dims,'r');
% end
% 
% 
% figure;
% hold on
% dims=[2 3];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% %plot(Zencl{4},dims);
% for iStep=1:timeSteps
%     plot(P{iStep},dims,'r');
% end
% 
% figure;
% hold on
% dims=[4 5];
% plot(Zencl{1},dims);
% plot(Zencl{2},dims);
% plot(Zencl{3},dims);
% %plot(Zencl{4},dims);
% for iStep=1:timeSteps
%     plot(P{iStep},dims,'r');
% end

%------------- END OF CODE --------------