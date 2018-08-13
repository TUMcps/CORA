function [nr] = select(obj,options,Rinit)
% select - selects the split strategy of the reachable set causing the
% least linearization error
%
% Syntax:  
%    [nr] = select(IHerrorActual,IHerrorAssume)
%
% Inputs:
%    IHerrorActual - cell array of actual linearizatuion errors
%    IHerrorAssume - assumed linearization error 
%
% Outputs:
%    nr - number of the selected split strategy
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      04-January-2008 
% Last update:  29-January-2008
%               29-June-2009
%               12-September-2017
% Last revision: ---

%------------- BEGIN CODE --------------

% compute all possible splits of the maximum reachable set
Rtmp=split(Rinit.set);
for i=1:length(Rtmp)
    Rxx{i}.set=Rtmp{i}{1}; %only test one of the two split sets
    %reset error
    Rxx{i}.error = 0*options.maxError;
end

if ~isempty(Rxx)

    % check performance index for all split sets
    perfInd = zeros(length(Rxx),1);
    for i=1:length(Rxx)
        % select reachable sets
        Rtest=Rxx{i};

        % obtain perf index
        [~,~,~,perfInd(i)] = linReach(obj,options,Rtest,0);
    end

    % find best performance index
    [~,ind]=sort(perfInd);
    nr=ind(1);
    
else
    %set to be split is empty
    nr = 0;
end



%------------- END OF CODE --------------