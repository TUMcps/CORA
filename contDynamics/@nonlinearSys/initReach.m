function [Rfirst, options] = initReach(obj,Rold,options)
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
%    Rfirst - first reachable set 
%    options - options struct
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
%               25-January-2016
%               11-August-2016
%               30-August-2017, AElguindy
%               12-September-2017, MA
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
    
    % set recursion
    if isfield(options,'recur')
        recur = options.recur;
    else
        recur = 1;
    end
    % compute reachable set of abstraction
    [Rti,Rtp,nr,perfInd] = linReach(obj,options,Rinit{iIteration},recur);


    %check if initial set has to be split
    if isempty(nr) %|| isempty(Rti)
        if recur == 1 || perfInd <= 1
            %save reachable sets in cell
            Rtotal_tp{setCounter}=Rtp;
            Rtotal_ti{setCounter}=Rti;
        else
            %save reachable sets in cell
            Rtotal_tp{setCounter}.set = zonotope([]);
            Rtotal_ti{setCounter} = zonotope([]);
        end
        
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
        
        
        Rres = initReach(obj,Rsplit,options);
        
        for i=1:length(Rres.tp)
            % add results to other results
            Rtotal_tp{setCounter} = Rres.tp{i};
            Rtotal_ti{setCounter} = Rres.ti{i};

            % update setCounter 
            setCounter = setCounter+1;
        end
    end
end


%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal_ti;

%------------- END OF CODE --------------