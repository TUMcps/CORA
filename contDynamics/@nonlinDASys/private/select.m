function dimForSplit = select(obj,options,Rinit,Rinit_y,iter)
% select - selects the split strategy of the reachable set causing the
%          least linearization error
%
% Syntax:
%    dimForSplit = select(obj,options,Rinit)
%
% Inputs:
%    obj - nonlinearSys or nonlinParamSys object
%    options - struct containing the algorithm settings
%    Rinit - initial reachable set (diff. variables)
%    Rinit_y - initial reachable set (alg. variables)
%    iter - flag for activating iteration
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

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       04-January-2008 
% Last update:   29-January-2008
%                29-June-2009
%                12-September-2017
%                02-January-2019 (NK, cleaned up the code)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute all possible splits of the maximum reachable set
Rtmp = split(Rinit);
R = cell(length(Rtmp),1);

for i = 1:length(Rtmp)
    R{i}.set = Rtmp{i}{1};        % only test one of the two split sets
    R{i}.error_x = zeros(size(options.maxError_x));
    R{i}.error_y = zeros(size(options.maxError_y));
end

% adapt the options for reachability analysis
options.maxError_x = Inf * options.maxError_x;

% loop over all split sets
perfInd = zeros(length(R),1);

for i=1:length(R)
    % compute the reachable set for the split set
    [~,~,~,perfInd(i)] = linReach(obj,options,R{i},Rinit_y,iter);
end

% find best performance index
[~,dimForSplit] = min(perfInd);

% ------------------------------ END OF CODE ------------------------------
