function [p]=plotReachProb(obj,options)
% plot - plots probabilistic reachable sets of a hybrid automaton
%
% Syntax:  
%    plotReach(obj)
%
% Inputs:
%    obj - hybrid automaton object
%
% Outputs:
%    p - probability of hitting unsafe set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 07-October-2007 
% Last update: 26-February-2008
% Last revision: ---

%------------- BEGIN CODE --------------


%compute projection matrix
dim=length(options.x0);
P=zeros(2,dim);
for iDim=1:2
    P(iDim,options.projectedDimensions(iDim))=1;
end

%load data from object structure
R=obj.result.reachSet.R;
loc=obj.result.reachSet.location;
dim=options.projectedDimensions;

%only CDC08prob!!-------------------------------------------------------------
%define unsafe set
IH=interval([-1e3; -1e3; -1e3; -1e3; -1e3], [1e3; -1.5; 1e3; 1e3; 1e3]);
IH=P*IH;
B=polytope(IH);
mArray1=[3,2.5,2,1.5,1,0.5]; %gammaSigma=2
mArray2=[2,1.7,1.5,1.2,1,0.8,0.5,0.2]; %gammaSigma=2
mArray3=linspace(options.gamma,0,10); %gammaSigma=2
%--------------------------------------------------------------------------

%plot each reachable set
for i=1:length(R.OT)
    for j=1:length(R.OT{i})
        %plot(Rred.p,dim); 
        if ~isempty(R.OT{i}{j})
            Rred=reduce(R.OT{i}{j},'girard',options.polytopeOrder);
            Rred=P*Rred;
            plot(Rred,'solid',[1 2],2.5); 
            %only CDC08prob------------------------------------------------
            %if i==2
                probTotal(j) = pyramid(Rred,mArray3,B);         
            %end
            %--------------------------------------------------------------
        end
        hold on
    end
end

%only CDC08prob
p=probTotal+1-erf(options.gamma/sqrt(2))^5;

%------------- END OF CODE --------------