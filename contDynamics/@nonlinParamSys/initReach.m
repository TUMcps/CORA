function [Rnext,options] = initReach(obj,Rinit,params,options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:
%    [Rnext,options] = initReach(obj,Rinit,params,options)
%
% Inputs:
%    obj - nonlinearSys object
%    Rinit - initial reachable set
%    params - model parameters
%    options - struct containing the algorithm settings
%
% Outputs:
%    Rfirst - first reachable set
%    options - struct containting the algorithm settings
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       29-October-2007 
% Last update:   04-January-2008
%                27-April-2009
%                16-August-2016 (harmonized with @nonlinearSys)
%                17-May-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization for the case that it is the first time step
    if ~iscell(Rinit)                                        
        R{1}.set = Rinit;
        R{1}.error = zeros(size(options.maxError));
        Rinit = R;
    end
    
    % compute reachable set using the options.linRem algorithm
    if strcmp(options.alg,'linRem')
        [Rnext,options] = aux_initReach_linRem(obj,Rinit,params,options);
        return;
    end

    % loop over all parallel sets
    setCounter = 1; Rtp = {}; Rti = {}; R0 = {};

    for i = 1:length(Rinit)

        % compute reachable set of abstraction
        [Rtemp_ti,Rtemp_tp,dimForSplit,opts] = linReach(obj,Rinit{i},params,options);
        
        % save prevErr to reuse
        if isfield(opts,'prevErr')
            options.prevErr = opts.prevErr;
        end
        
        % check if initial set has to be split
        if isempty(dimForSplit)
            Rtp{setCounter} = Rtemp_tp;
            Rtp{setCounter}.prev = i;
            Rti{setCounter} = Rtemp_ti;
            R0{setCounter} = Rinit{i};
            setCounter = setCounter+1;
        else
            disp('split!!');

            % split the initial set 
            Rtmp = split(Rinit{i}.set,dimForSplit);
            Rsplit{1}.set = Rtmp{1};
            Rsplit{2}.set = Rtmp{2};

            % reset the linearization error
            Rsplit{1}.error = 0*options.maxError;
            Rsplit{2}.error = 0*options.maxError;

            % compute the reachable set for the splitted sets
            Rres = initReach(obj,Rsplit,params,options);

            for j = 1:length(Rres.tp)
                Rtp{setCounter} = Rres.tp{j};
                Rtp{setCounter}.parent = i;
                Rti{setCounter} = Rres.ti{j};
                R0{setCounter} = Rres.R0{j};
                setCounter = setCounter + 1;
            end
        end
    end

    % store the results
    Rnext.tp = Rtp;
    Rnext.ti = Rti;
    Rnext.R0 = R0;
end


% Auxiliary functions -----------------------------------------------------

function [Rnext,options] = aux_initReach_linRem(obj,Rinit,params,options)

    % compute the reachable set using the options.alg = 'lin' algorithm to
    % obtain a first rough over-approximation of the reachable set
    options_ = options;
    options_.alg = 'lin';
    
    Rtemp = initReach(obj,Rinit,params,options_);
    
    % loop over all parallel sets 
    Rtp = cell(length(Rtemp.ti),1); 
    Rti = cell(length(Rtemp.ti),1);

    for i = 1:length(Rtemp.ti)

        % compute reachable set with the "linear remainder" algorithm to
        % obtain a refined reachable set
        R0 = Rtemp.R0{i};
        options.Ronestep = Rtemp.ti{i};
        
        [Rti{i},Rtp{i}] = linReach(obj,R0,params,options);
        
        % copy information about previous reachble set and parent
        Rtp{i}.prev = Rtemp.tp{i}.prev;
        if isfield(Rtemp.tp{i},'parent')
           Rtp{i}.parent = Rtemp.tp{i}.parent; 
        end
    end

    % store the results
    Rnext.tp = Rtp;
    Rnext.ti = Rti;
end

% ------------------------------ END OF CODE ------------------------------
