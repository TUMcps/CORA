function [Rfirst,options] = initReach(obj,Rold,options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:  
%    [Rnext] = initReach(obj,R,options)
%
% Inputs:
%    obj - nonlinearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - nonlinearSys object
%    Rfirst - first reachable set 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      29-October-2007 
% Last update:  04-January-2008
%               27-April-2009
%               16-August-2016 (now identical to initReach of @nonlinearSys; solve by inheritance in the future)
% Last revision: ---

%------------- BEGIN CODE --------------

%check if Rinit can be obtained from options
if iscell(Rold)
    Rinit=Rold;
else
    %first run
    Rinit{1}.set=Rold;
    Rinit{1}.error=0*options.maxError; %initially error is assumed to be 0
end
iterations=length(Rinit);

%initialize set counter
setCounter=1;

for iIteration=1:iterations
    %try
    [Rti,Rtp,nr,options] = linReach(obj,options,Rinit{iIteration},1);
%     catch
%         disp('error')
%     end

    %check if initial set has to be split
    if isempty(nr)
        %save reachable sets in cell
        Rtotal_tp{setCounter}=Rtp;
        Rtotal_ti{setCounter}=Rti;
        
        %setCounter update
        setCounter=setCounter+1;
    else
        disp('split!!');

        %split initial set 
        Rtmp=split(Rinit{iIteration}.set,nr);
        Rsplit{1}.set=Rtmp{1};
        Rsplit{2}.set=Rtmp{2};
        %reset error
        Rsplit{1}.error = 0*options.maxError;
        Rsplit{2}.error = 0*options.maxError;
        
%         figure
%         plot(Rinit{iIteration});
%         plot(Rsplit{1});
%         plot(Rsplit{2});
        
        [Rres,options]=initReach(obj,Rsplit,options);
        
        for i=1:length(Rres.tp)
            % add results to other results
            Rtotal_tp{setCounter}=Rres.tp{i};
            Rtotal_ti{setCounter}=Rres.ti{i};

            % update setCounter 
            setCounter=setCounter+1;
        end
    end
end


%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal_ti;

%------------- END OF CODE --------------