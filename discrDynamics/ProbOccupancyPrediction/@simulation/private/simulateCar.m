function [pA,pTotal]=simulateCar(runs,T,Omega,pA0,pB,Tx,pX,RoW,led,qFree,xSegment)
% simulateCar - simulates a single car, based on the transition matrix T,
% the interaction matrix Omega, the probability pB of the vehicle ahead, 
% the transition matrix Tx for the change from approaching to crossing, 
% the RoW (right of way) flag, and the segment where the RoW path is 
% intersected (xSegment) 
%
% Syntax:  
%    [pA,pX]=simulateCar(T,Omega,pB,Tx,pX,RoW,xSegment)
%
% Inputs:
%    runs - number of computed time intervals
%    T - transition matrix
%    Omega - interaction matrix
%    pB - probability vector of the car ahead
%    Tx - transition matrix from approaching to crossing
%    pX - probability that a car is crossing that has RoW
%    RoW - 1/0; 1: car has right of way; 0: otherwise
%    led - 1/0; 1: car is leading vehicle; 0: otherwise
%    qFree - probability distribution of driv
%    xSegment - segment, where path is crossed
%
% Outputs:
%    pA - probability distribution of car A
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 16-November-2007 
% Last update: 13-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

%determine number of acceleration modes
modes=length(Omega);
if modes==0
    modes=4;
end


%define zero vector z and identity matrix I
z=0*pA0{1}{1}{1}.T;
I=eye(length(z));

%initialize pA1total, pA2total for first time step/interval
pA1total{1}.T=z; 
% it is assumed, that probability of beeing in crossing mode is zero at the 
% beginning
pA2total{1}.T=z;

for iMode=1:modes
    pA1total{1}.T=pA1total{1}.T+pA0{1}{1}{iMode}.T;
end
pA1total{1}.OT=pA1total{1}.T;
pA2total{1}.OT=pA2total{1}.T;

%complement pX
if ~RoW
    for i=1:runs
%         newpX{i}=z;
%         newpX{i}(xSegment+1)=pX(i);
        newpXappr{i}=z;
        newpXappr{i}(xSegment+1)=1;
    end
%     pX=newpX;
    pXappr=newpXappr; %virtual car of probability 1 for approaching mode
end


%set first value for pA
pA=pA0;

for iStep=1:runs
    
    % compute mode matrix if it is no leading vehicle
    if ~(led & RoW)
        % initialize total mode vector
        M1total=z;
        M2total=z;

        for iMode=1:modes
            %initialize M vector
            M1{iMode}=z;
            M1x{iMode}=z;
            M2{iMode}=z;
            M2x{iMode}=z;

            for iAheadMode=1:modes
                if RoW %<-- continue here
                    %compute unnormalized M vector
                    M1{iMode}=M1{iMode}+Omega{iMode,iAheadMode}*pB{iStep}{1}{iAheadMode}.T; %A drives after B             
                else
                    if led
                        %compute unnormalized M vector
                        M1{iMode}=qFree(iMode)*ones(length(z),1); %<-- changed
                        %crossing mode
                        M2{iMode}=qFree(iMode)*ones(length(z),1); %<-- changed                          
                    else
                        %compute unnormalized M vector
                        M1{iMode}=M1{iMode}+Omega{iMode,iAheadMode}*pB{iStep}{1}{iAheadMode}.T; %A drives after B 
                        %crossing mode
                        M2{iMode}=M2{iMode}+Omega{iMode,iAheadMode}*pB{iStep}{2}{iAheadMode}.T; %A drives in crossing mode after B
                    end
                end
            end
            if ~RoW
                M1x{iMode}=M1x{iMode}+Omega{iMode,1}*pXappr{iStep}; %A drives after X
                M1{iMode}=M1{iMode}.*M1x{iMode};   
            end
            %correct due to probabilities that represent that the car has been
            %overtaken
            M1total=M1total+M1{iMode};    
            M2total=M2total+M2{iMode};  
        end


        M1total=M1total+1e-4; %avoid division by 0
        M2total=M2total+1e-4; %avoid division by 0
        for iMode=1:modes
            M1{iMode}=M1{iMode}./M1total; %normalize
            M2{iMode}=M2{iMode}./M2total; %normalize
        end
    end
    
    %initialize pA1total, pA2total
    pA1total{iStep+1}.T=z;
    pA2total{iStep+1}.T=z;
 
    pA1total{iStep}.OT=z; %<--changed from iStep+1 to iStep!!
    pA2total{iStep}.OT=z; %<--changed from iStep+1 to iStep!!
    
    for iMode=1:modes
        %if car is leading vehicle
        if led & RoW
            if mod(iStep,1)==0 %<-- change updatefrequency here
                pA{iStep+1}{1}{iMode}.T=T.T{iMode}*qFree(iMode)*pA1total{iStep}.T;
                pA{iStep}{1}{iMode}.OT=T.OT{iMode}*qFree(iMode)*pA1total{iStep}.T;   %<--changed from iStep+1 to iStep!!          
            else
                pA{iStep+1}{1}{iMode}.T=T.T{iMode}*pA{iStep}{1}{iMode}.T;
                pA{iStep}{1}{iMode}.OT=T.OT{iMode}*pA{iStep}{1}{iMode}.T; %<--changed from iStep+1 to iStep!!
            end
            
            %normalize (.T values is enough)
            %pA{iStep+1}{1}{iMode}.T=normalize(pA{iStep+1}{1}{iMode}.T,...
            %    pA{iStep}{1}{iMode}.T);   %<--normalization commented
            
        %if car is not leading vehicle and has right of way
        elseif RoW
            pTemp=diag(M1{iMode})*pA1total{iStep}.T;
            pA{iStep+1}{1}{iMode}.T=T.T{iMode}*pTemp; 
            pA{iStep}{1}{iMode}.OT=T.OT{iMode}*pTemp; %<--changed from iStep+1 to iStep!!
            
            %normalize (.T values is enough)
            %pA{iStep+1}{1}{iMode}.T=normalize(pA{iStep+1}{1}{iMode}.T,pTemp); %<--normalization commented
            
        %if car has not right of way     
        else
            pTemp1=diag(M1{iMode})*pA1total{iStep}.T;
            pTemp2=diag(M2{iMode})*pA2total{iStep}.T;

            
            %phase 1: approaching
            %time point
            pA{iStep+1}{1}{iMode}.T=T.T{iMode}*pTemp1;
            %normalize (.T values is enough)
            %pA{iStep+1}{1}{iMode}.T=normalize(pA{iStep+1}{1}{iMode}.T,pTemp1); 
            %from approaching to crossing
            pAnew=(I-(1-pX(iStep))*Tx)*pA{iStep+1}{1}{iMode}.T;
            
            %phase 2: crossing
            pA{iStep+1}{2}{iMode}.T=T.T{iMode}*pTemp2;
            %normalize (.T values is enough)
            %pA{iStep+1}{2}{iMode}.T=normalize(pA{iStep+1}{2}{iMode}.T,pTemp2);  
            %from approaching to crossing
            pA{iStep+1}{2}{iMode}.T=pA{iStep+1}{2}{iMode}.T+(1-pX(iStep))*Tx*pA{iStep+1}{1}{iMode}.T;
            
            %Important change: update approaching after crossing
            %computations
            pA{iStep+1}{1}{iMode}.T=pAnew;
            
            %time interval
            pA{iStep}{1}{iMode}.OT=T.OT{iMode}*pTemp1; %<--changed from iStep+1 to iStep!!
            pA{iStep}{1}{iMode}.OT=(I-(1-pX(iStep))*Tx)*pA{iStep}{1}{iMode}.OT; %<--changed from iStep+1 to iStep!!
            pA{iStep}{2}{iMode}.OT=T.OT{iMode}*pTemp2; %<--changed from iStep+1 to iStep!!
            pA{iStep}{2}{iMode}.OT=pA{iStep}{2}{iMode}.OT+(1-pX(iStep))*Tx*pA{iStep}{1}{iMode}.OT;  %<--changed from iStep+1 to iStep!!                                     
        end
        

        %compute total probabilities 
        pA1total{iStep+1}.T=pA1total{iStep+1}.T+pA{iStep+1}{1}{iMode}.T;
        pA1total{iStep}.OT=pA1total{iStep}.OT+pA{iStep}{1}{iMode}.OT; %<--changed from iStep+1 to iStep!!
        if ~RoW
            %compute total probabilities for crossing mode
            pA2total{iStep+1}.T=pA2total{iStep+1}.T+pA{iStep+1}{2}{iMode}.T;
            pA2total{iStep}.OT=pA2total{iStep}.OT+pA{iStep}{2}{iMode}.OT; %<--changed from iStep+1 to iStep!!
        end
       
    end
    %define pTotal
    pTotal.T{iStep}=pA1total{iStep}.T;   
    pTotal.OT{iStep}=pA1total{iStep}.OT;
    if ~RoW
        %add total probabilities from crossing mode
        pTotal.T{iStep}=pTotal.T{iStep}+pA2total{iStep}.T;   
        pTotal.OT{iStep}=pTotal.OT{iStep}+pA2total{iStep}.OT;
    end
end


%------------- END OF CODE --------------