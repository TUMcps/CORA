function [crash] = crashCheckOld(carAstate,carBstate,carAacc,carBacc,T)
% crashCheck - checks if car A can crash into car B if both cars accelerate
% with the given accelerations carAacc and carBacc for a time horizin T
% and an emergency brake with maxAcc of both cars afterwards.
%
% Syntax:  
%    [crash] = (carAstate,carBstate,carAacc,carBacc,T,maxAcc)
%
% Inputs:
%    carAstate - initial state of car A
%    carBstate - initial state of car B
%    carAacc - acceleration of car A
%    carBacc - acceleration of car B
%    T - time interval for carAacc, carBacc 
%
% Outputs:
%    crash - 1 if crash occurs and 0 otherwise
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 26-October-2007 
% Last update: 27-June-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%set maximum acceleration
maxAcc=-10;

%compute stopping times 
stopTimeA=getStopTime(carAstate,carAacc,T,maxAcc);
stopTimeB=getStopTime(carBstate,carBacc,T,maxAcc);

%sort stopping times and step size T
times=sort([0,stopTimeA,stopTimeB,T]);

%get states of sorted times
[xValueA,vValueA,accValueA]=getStates(carAstate,carAacc,T,maxAcc,times);
[xValueB,vValueB,accValueB]=getStates(carBstate,carBacc,T,maxAcc,times);

%check for crash in all time intervals
crash=check(xValueA,vValueA,accValueA,...
    xValueB,vValueB,accValueB,times);
end



%check if a crsh occurs---------------------------------------
function crash=check(xValueA,vValueA,accValueA,...
    xValueB,vValueB,accValueB,times)
    
    %init crash
    crash=0;    
    
    %loop
    for iTime=1:(length(xValueA)-1)
        
        %set values
        xA=xValueA(iTime);
        vA=vValueA(iTime);
        aA=accValueA(iTime);
        xB=xValueB(iTime);
        vB=vValueB(iTime);
        aB=accValueB(iTime);
        
        %compute determinant
        determinant=(vB-vA)^2-2*(aB-aA)*(xB-xA);
    
        %check if determinant is positive
        if determinant>=0
            
            %obtain delta time
            deltaTime=times(iTime+1)-times(iTime);
            
            if (aA==aB) & (abs(vA-vB)<1e-6)
                %no crash
            elseif (aA==aB) 
                tCrash=(xA-xB)/vB-vA;

                %crash time in time interval?
                if (tCrash>0) & (tCrash<deltaTime)
                    crash=1;
                end
            else
                tCrash1=((vA-vB)+sqrt(determinant))/(2*(aB-aA));
                tCrash2=((vA-vB)-sqrt(determinant))/(2*(aB-aA));

                %crash time in time interval?
                if (tCrash1>0) & (tCrash1<deltaTime)
                    crash=1;
                elseif (tCrash2>0) & (tCrash2<deltaTime)
                    crash=1;                  
                end            
            end
        end
    end
end


%simulate a single vehicle until it stops---------------------
function stopTime=getStopTime(state,acc,T,maxAcc)

    %assumed stopping time for acc
    assumedStop=-state(2)/acc;

    if (assumedStop>=0) & (T>=assumedStop)
        %stop time equals assumed stop
        stopTime=assumedStop;
    else
        %position and velocity after time T
        xT=0.5*acc*T^2+state(2)*T+state(1);
        vT=acc*T+state(2);

        %compute stop time
        stopTime=T-vT/maxAcc;
    end
end

%simulate a single vehicle until it stops---------------------
function [xValue,vValue,accValue]=getStates(state,acc,T,maxAcc,times)

    %init values
    xValue(1)=state(1);
    vValue(1)=state(2);
    accValue(1)=acc;
    
    for iTime=2:length(times)
        %get time
        t=times(iTime)-times(iTime-1);

        %compute next position and velocity
        xValue(iTime)=0.5*accValue(iTime-1)*t^2+vValue(iTime-1)*t+xValue(iTime-1);
        vValue(iTime)=accValue(iTime-1)*t+vValue(iTime-1);
        
        %compute acceleration values
        if vValue(iTime)>0
            %compute acceleration value if t<T
            if times(iTime)<T                
                accValue(iTime)=acc;
            %compute acceleration values if t>=T
            else
                accValue(iTime)=maxAcc;
            end
        else
            accValue(iTime)=0;
        end
    end
end


%------------- END OF CODE --------------