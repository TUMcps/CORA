function pred = calibrate(pred,data,n_out,options)
% calibrate - calibrate uncertainties of set-based predictor
%
% Syntax:
%    pred = calibrate(pred,data,options)
%    pred = calibrate(pred,data,options,n_out)
%
% Inputs:
%    pred - confPredictor object
%    data - calibration data, either specified as trajectory objects (where
%           traj(i).x is the i'th input and traj(i).y the corresponding 
%           output) or as struct with data.input(:,i) and data.target(:,i)
%    n_out - number of data points allowed outside of prediction set
%    options - specifications
%
% Outputs:
%    pred - confPredictor object with calibrated uncertainty set
%
% References: 
%    [1] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based 
%        Uncertainty Quantification for Regression and Classification 
%        Tasks," arXiv, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: confPredictor

% Authors:       Laura Luetzow
% Written:       25-September-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 3
    n_out = 0;
end
if nargin < 4
    options = struct();
end

% add data to the predictor object
if ~isa(data,'trajectory')
    % transform to trajectory objects
    traj = [];
    u = zeros(pred.sys.nrOfInputs,1);
    for i = 1:size(data.input,2)
        traj = [traj; trajectory(u,data.input(:,i),data.target(:,i),[],1)];
    end
    pred.testSuite = traj;
else
    pred.testSuite = data;
end

% calibrate uncertainties
if strcmp(pred.type,'conformalAdd')
    % compute additive uncertainty for each dimension with conformal
    % prediction
    pred = aux_calibrateConformalAdd(pred, n_out);
else
    % select uncertainties from all candidate uncertainties and calibrate
    pred = aux_calibrateOthers(pred, n_out, options);
end
end


% Auxiliary functions -----------------------------------------------------

function pred = aux_calibrateConformalAdd(pred, n_out)
% compute additive uncertainty for each dimension with conformal
% prediction

% compute quantile
num_samples = length(pred.testSuite);
p_quantile = n_out / num_samples;

if strcmp(pred.task,'class')
    % classification tasks
    s = zeros(1,num_samples);
    for m=1:length(pred.testSuite)
        % sort softmax predictions in descending order
        y_pred = pred.sys.out_mFile(pred.testSuite(m).x, pred.testSuite(m).u);
        y_softmax = softmax(y_pred);
        y_true = pred.testSuite(m).y;
        s(:,m) = 1 - y_softmax(logical(y_true));
    end
else
    % regression tasks
    s = zeros(pred.sys.nrOfOutputs,num_samples);
    for m=1:length(pred.testSuite)
        y_pred = pred.sys.out_mFile(pred.testSuite(m).x, pred.testSuite(m).u);
        s(:,m) = abs(pred.testSuite(m).y - y_pred);
    end
end
q = quantile(s,1-p_quantile,2);

% establish conformance
if strcmp(pred.task,'class')
    pred.U = q;
else
    pred.U = interval(-q,q);
end
pred.results.fval = sum(2*q);
pred.results.p = [q];
end

function pred = aux_calibrateOthers(pred,n_out,options)
if isempty(pred.U_init)
    if isfield(options,'up_spec')
        pred = initializeUncertainties(pred,options.up_spec);
    else % use default (identity matrix)
        pred = initializeUncertainties(pred);
    end
end

% construct parameter objects for conformance synthesis
if ~isa(pred.U_init, 'zonotope')
    % convert uncertainty set to a zonotope
    U_init = zonotope(pred.U_init);
else
    U_init = pred.U_init;
end
params.U = U_init;
params.R0 = zonotope(zeros(pred.nrOfInputs,1));
params.testSuite = pred.testSuite;
params.tFinal = 0;

% set options
options.cs.task = pred.task;
options.cs.numOutlier = n_out;

% establish conformance
% outlier detection
[params, results] = conform(pred.sys,params,options);

if strcmp(pred.type,'interval')
    pred.U = interval(params.U);
else
    pred.U = params.U;
end
pred.results = results;
end

% ------------------------------ END OF CODE ------------------------------
