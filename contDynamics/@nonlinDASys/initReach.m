function [Rnext, options] = initReach(nlnsysDA, Rstart, params, options)
% initReach - computes the reachable continuous set for the next step
%    given the current set as the start set
%
% Syntax:
%    [Rnext, options] = initReach(nlnsysDA, Rstart, params, options)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%    Rstart - initial reachable set
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - nonlinearSys object
%    Rnext - next reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       22-November-2011
% Last update:   08-August-2016
%                11-January-2021 (MW, enable splitting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% regular iteration (incl. splits)
if isstruct(Rstart)
    Rinit = Rstart.tp;
    Rinit_y = Rstart.y;
    
    iterations=length(Rinit);
    
% initialization for first time step / after split
else
    
    % nr of init sets can vary if split in first step
    if ~iscell(Rstart)
        Rstart = {Rstart};
    end
    iterations = length(Rstart);
    
    Rinit = cell(iterations,1); Rinit_y = cell(iterations,1);
    for k=1:iterations
        Rinit{k}.set = Rstart{k};
        Rinit{k}.error_x = zeros(size(options.maxError_x));
        Rinit{k}.error_y = zeros(size(options.maxError_y));
        % obtain consistent initial algebraic set
        y0 = params.y0guess;
        y0 = consistentInitialState(nlnsysDA, center(Rinit{k}.set), y0, params.uTrans);
        Rinit_y{k} = zonotope(y0);
    end

end

%initialize set counter
setCounter=1;

for k=1:iterations

    [Rti,Rtp,Rti_y,~,dimForSplit,options] = linReach(nlnsysDA,Rinit{k},Rinit_y{k},params,options);

    %check if initial set has to be split
    if isempty(dimForSplit)
        %save reachable sets in cell
        Rtotal_tp{setCounter} = Rtp;
        Rtotal_tp{setCounter}.prev = k;
        Rtotal_ti{setCounter} = Rti;
        Rtotal_y{setCounter} = Rti_y;
        %setCounter update
        setCounter = setCounter+1;
    else
        disp('split!!');

        % split initial set 
        Rsplit = split(Rinit{k}.set,dimForSplit);
        % adapt y0guess given current Rinit_y
        if options.t > options.tStart
            options.y0guess = center(Rinit_y{1});
        end
        
        [Rres,options] = initReach(nlnsysDA,Rsplit,params,options);
        
        for i=1:length(Rres.tp)
            % add results to other results
            Rtotal_tp{setCounter} = Rres.tp{i};
            Rtotal_tp{setCounter}.parent = k;
            Rtotal_ti{setCounter} = Rres.ti{i};
            Rtotal_y{setCounter} = Rres.y{i};
            % update setCounter 
            setCounter=setCounter+1;
        end
    end
end

% write results to reachable set struct Rfirst
Rnext.tp = Rtotal_tp;
Rnext.ti = Rtotal_ti;
Rnext.y = Rtotal_y;

% ------------------------------ END OF CODE ------------------------------
