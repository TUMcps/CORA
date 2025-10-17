function Y = predict(pred,input)
% predict - contruct prediction set for the given input
%
% Syntax:
%    Y = predict(pred,input)
%
% Inputs:
%    pred - confPredictor object
%    input - dim_x x 1 input vector
%
% Outputs:
%    Y- prediction set
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

% compute the prediction set using linearization
u_ref = zeros(pred.sys.nrOfInputs,1);
f0 = pred.sys.out_mFile(input, u_ref); % point prediction

if strcmp(pred.type,'conformalAdd')
    % add uncertainty for each dimension computed with conformal
    % prediction
    if isa(pred.U,'interval')
        % for regression tasks
        Y = f0 + pred.U;
    else
        % for classification tasks
        y = softmax(f0);
        idz = 1:length(y);
        Y = idz(y >= 1-pred.U-1e-15); % subtract 1e-15 for numerical robustness
    end
else
    % construct prediction set using the Jacobian matrix
    [~,D_bar] = pred.sys.out_jacobian(input, u_ref);
    Y = f0 + D_bar*pred.U; % prediction set
end

if strcmp(pred.type,'interval')
    % overapproximate prediction set with an interval
    Y = interval(Y);
end
end

% ------------------------------ END OF CODE ------------------------------
