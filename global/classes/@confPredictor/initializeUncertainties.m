function pred = initializeUncertainties(pred,up_spec)
% initializeUncertainties - set the initial uncertainty set uncertainties 
%   of a set-based predictor
%
% Syntax:
%    pred = initializeUncertainties(pred)
%    pred = initializeUncertainties(pred,up_spec)
%
% Inputs:
%    pred - confPredictor object
%    up_spec - specifications for uncertainty placement
%       .G_u - generator template for the uncertainty set (if not
%           specified: set to identity matrix)
%       .p_nonzero - boolean vector specifying which uncertainties will be
%           identified, while the remaining ones are set to zero 
%       .method - placement strategy to determine p_nonzero if not
%           specified
%       .numUnc - number of uncertainties to be identified
%       .addRandGens - boolean determining if additional random generators
%           will be added
%
% Outputs:
%    pred - confPredictor object with updated initial uncertainty set
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
% Written:       29-September-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(pred.type, "conformalAdd")
    throw(CORAerror("CORA:specialError", "Predictor of type " + ...
        "conformalAdd does not require initialization of uncertainty set."))
end

% construct uncertainty set template
if nargin > 1 && isfield(up_spec,'G_u')
    G_u = up_spec.G_u;
else
    % set generator matrix to identity matrix
    G_u = eye(pred.nrOfUncertainties);
end

if nargin > 1
    % place uncertainties according to up_spec if specified
    if isfield(up_spec,'p_nonzero')
        p_nonzero = up_spec.p_nonzero;
    elseif isfield(up_spec,'method') && isfield(up_spec,'numUnc')
        numUnc = up_spec.numUnc;
        % up_spec describes uncertainty placement method -> use
        % method to determine the uncertainties to identify
        p_nonzero = aux_placeUncertainties(pred,up_spec.method,numUnc);
    else
        % identify all uncertainties
        p_nonzero = true(size(G_u,1),1);
    end

    if isfield(up_spec,'addRandGens') && up_spec.addRandGens
        % add random generators to G_u
        G_u = [G_u randn(size(G_u,1), sum(p_nonzero))];
    end

    % adapt G_u to only have generators pointing in the
    % directions of the uncertainties to identify
    G_u(~logical(p_nonzero),:) = 0;
    G_u(:,sum(abs(G_u),1)==0) = [];
end
U_init = zonotope(zeros(pred.nrOfUncertainties,1), G_u);

if strcmp(pred.type, "zonotope")
    pred.U_init = U_init;
elseif strcmp(pred.type, "interval")
    pred.U_init = interval(U_init);
else
    throw(CORAerror("CORA:specialError", "Initialization not " + ...
        "implemented for predictor type."))
end
end


% Auxiliary functions -----------------------------------------------------


function p_nonzero = aux_placeUncertainties(pred,method,numUnc)
% select uncertainties that should be identify

% number of parametric uncertainties
n_p = pred.nrOfUncertainties - pred.nrOfOutputs;
n_y = pred.nrOfOutputs;

% total number of identified uncertainties
n_p_id = numUnc - n_y;

% determine where to place the uncertainties
switch method
    case "all"
        % identify all uncertainties
        p_nonzero = ones(pred.sys.nrOfInputs,1);

        if n_p_id < n_p
            p_nonzero(1:end-numUnc) = 0;
        end

    case "onlyOutput"
        % only output uncertainties
        p_nonzero = zeros(pred.sys.nrOfInputs,1);
        p_nonzero(n_p+1:end) = 1;

    case "onlyInput"
        % only input uncertainties
        p_nonzero = zeros(pred.sys.nrOfInputs,1);
        p_nonzero(1:pred.sys.nrOfDims) = 1;

    case "random"
        % random uncertainties

        % identify n_p_id random parametric uncertainties
        p_nonzero = zeros(pred.sys.nrOfInputs,1);
        helper = randperm(n_p+n_y, n_p+n_y);
        helper = helper(1:n_p_id+n_y);
        p_nonzero(helper) = 1;

    case "outRandom"
        % random uncertainties

        % % identify all output uncertainties
        p_nonzero = zeros(pred.sys.nrOfInputs,1);
        p_nonzero(n_p+1:end) = 1;

        % identify n_p_id random parametric uncertainties
        helper = randperm(n_p, n_p);
        helper = helper(1:n_p_id);
        p_nonzero(helper) = 1;

    case "QR"
        % % identify all output uncertainties
        p_nonzero = zeros(pred.sys.nrOfInputs,1);
        p_nonzero(n_p+1:end) = 1;

        % identify n_p_id parametric uncertainties selected using QR factorization
        D_2d = zeros(pred.sys.nrOfOutputs*length(pred.testSuite),pred.sys.nrOfInputs);
        for i_d = 1: length(pred.testSuite)
            % compute prediction error and sensitivty
            [~, D_2d((i_d-1)*pred.sys.nrOfOutputs+1:i_d*pred.sys.nrOfOutputs,:)] = pred.sys.out_jacobian(pred.testSuite(i_d).x, zeros(pred.sys.nrOfInputs,1));
        end
        [~,~,P] = qr(D_2d,"vector");
        P(P>n_p) = [];
        p_nonzero(P(1:n_p_id)) = 1;
end
end

% ------------------------------ END OF CODE ------------------------------
