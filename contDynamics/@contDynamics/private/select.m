function dimForSplit = select(obj,options,Rinit)
% select - selects the split strategy of the reachable set causing the
%          least linearization error
%
% Syntax:  
%    dimForSplit = select(obj,options,Rinit)
%
% Inputs:
%    obj - nonlinearSys or nonlinParamSys object
%    options - struct containing the algorithm settings
%    Rinit - initial reachable set
%
% Outputs:
%    dimForSplit - dimension that is split to reduce the linearization
%                  error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linReach

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      04-January-2008 
% Last update:  29-January-2008
%               29-June-2009
%               12-September-2017
%               02-January-2019 (NK, cleaned up the code)
% Last revision: ---

%------------- BEGIN CODE --------------

% compute all possible splits of the maximum reachable set
Rtmp = split(Rinit.set);
R = cell(length(Rtmp),1);

for i = 1:length(Rtmp)
    R{i}.set = Rtmp{i}{1};        % only test one of the two split sets
    R{i}.error = zeros(size(options.maxError));
end

% adapt the options for reachability analysis
maxError = options.maxError;
options.maxError = inf * maxError;

if strcmp(options.alg,'linRem')
   options.alg = 'lin'; 
end

% loop over all split sets
perfInd = zeros(length(R),1);

for i=1:length(R)
    
    % compute the reachable set for the splitted set
    [~,Rtp] = linReach(obj,options,R{i});
    
    % compute performance index (max lin. error) for the split 
    perfInd(i) = max(Rtp.error./maxError);
end

% find best performance index
[~,dimForSplit] = min(perfInd);

%------------- END OF CODE --------------