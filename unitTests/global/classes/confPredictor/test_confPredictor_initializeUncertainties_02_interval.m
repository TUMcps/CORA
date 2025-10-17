function res = test_confPredictor_initializeUncertainties_02_interval
% test_confPredictor_initializeUncertainties_02_interval - unit test for
%       uncertainty placement for interval-based predictors
%
% Syntax:
%    res = test_confPredictor_initializeUncertainties_02_interval
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% user specifications

% define predictor model with many candidate uncertainties
regr_func = @(x) [ones(size(x)); x; x.^2; x.^3; x.^4];
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(regr_func(x).*u(1:10,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(regr_func(x).*u(11:20,:),1)];
dim_x = 2;
dim_u = 20;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create zono-conformal predictor 
pred = confPredictor(sys,"interval");

% use default generator template
G_u_def = eye(dim_u); % default is identity matrix
pred = initializeUncertainties(pred);
assert(isa(pred.U_init,'interval'))
Z = zonotope(pred.U_init);
assert(all(Z.G == G_u_def,'all'))

% use predefined generator template
up_spec1.G_u = randn(dim_u,15);
pred = initializeUncertainties(pred,up_spec1);
assert(isa(pred.U_init,'interval'))
assert(all(pred.U_init.sup-pred.U_init.inf> 1e-6)) 

% predefine the nonzero uncertainties
up_spec2.p_nonzero = [1 1 zeros(1,dim_u-2)]';
pred = initializeUncertainties(pred,up_spec2);
assert(isa(pred.U_init,'interval'))
diffs = pred.U_init.sup-pred.U_init.inf;
assert(all(diffs(logical(up_spec2.p_nonzero)) > 1e-6)) % nonzero dimensions
assert(all(diffs(~logical(up_spec2.p_nonzero)) < 1e-6)) % zero dimensions

% test different uncertainty placement strategies with different numbers of
% uncertainties
pred_all = initializeUncertainties(pred);
for numUnc = [5 10 15 20]
    up_spec3.numUnc = numUnc;
    for up_method = ["random","QR","outRandom","outRandom*"]
        if strcmp(up_method,"outRandom*")
            % add random generator to generator matrix G_U
            up_spec3.method = "outRandom";
            up_spec3.addRandGens = true;
        else
            up_spec3.method = up_method;
            up_spec3.addRandGens = false;
        end
        pred = initializeUncertainties(pred,up_spec3);
        assert(isa(pred.U_init,'interval'))

        if numUnc == dim_u && ~strcmp(up_method,"outRandom*")
            % check if same results for identifying all uncertainties
            assert(all(withinTol(pred_all.U_init.sup,pred.U_init.sup,1e-6),'all'))
            assert(all(withinTol(pred_all.U_init.inf,pred.U_init.inf,1e-6),'all'))
        end
    end
end

res = true;
end

% ------------------------------ END OF CODE ------------------------------
