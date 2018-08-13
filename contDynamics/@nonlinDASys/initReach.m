function [Rfirst, options] = initReach(obj, Rold, options)
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
% Written:      22-November-2011
% Last update:  08-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%check if Rinit can be obtained from options
if isfield(Rold,'tp') && iscell(Rold.tp)
    Rinit = Rold.tp;
    Rinit_y = Rold.y;
else
    Rinit{1} = Rold;
    % obtain consistent initial algebraic set
    y0 = options.y0guess;
    y0 = consistentInitialState(obj, options.x0, y0, options.uTrans);
    Rinit_y{1} = zonotope(y0);
    %initialize old error (should be removed when adapting concept of @nonlinearSys)
    options.oldError = 0*options.maxError;
    options.oldError_x = 0*options.maxError_x;
    options.oldError_y = 0*options.maxError_y;
end
iterations=length(Rinit);

%initialize set counter
setCounter=1;

for iIteration=1:iterations

    [Rti,Rtp,Rti_y,perfInd,nr,options] = linReach(obj,options,Rinit{iIteration},Rinit_y{iIteration},1);

    %check if initial set has to be split
    if isempty(nr)
        %save reachable sets in cell
        Rtotal_tp{setCounter} = Rtp;
        Rtotal_ti{setCounter} = Rti;
        Rtotal_y{setCounter} = Rti_y;
        %setCounter update
        setCounter=setCounter+1;
    else
        disp('split!!');

        %split initial set 
        Rsplit=split(Rinit{iIteration},options,nr);
        
%         figure
%         plot(Rinit{iIteration});
%         plot(Rsplit{1});
%         plot(Rsplit{2});
        
        [obj,Rres,options]=initReach(obj,Rsplit,options);
        
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
Rfirst.y = Rtotal_y;

%------------- END OF CODE --------------