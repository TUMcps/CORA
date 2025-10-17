function res = test_confPredictor_initializeUncertainties_01_zcp
% test_confPredictor_initializeUncertainties_01_zcp - unit test for
%       uncertainty placement for zono-conformal predictors
%
% Syntax:
%    res = test_confPredictor_initializeUncertainties_01_zcp
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
pred = confPredictor(sys,"zonotope");

% use default generator template
G_u_def = eye(dim_u); % default is identity matrix
pred = initializeUncertainties(pred);
assert(all(pred.U_init.G == G_u_def,'all'))
assert(isa(pred.U_init,'zonotope'))

% use predefined generator template
up_spec1.G_u = randn(dim_u,15);
pred = initializeUncertainties(pred,up_spec1);
assert(all(pred.U_init.G == up_spec1.G_u,'all'))
assert(isa(pred.U_init,'zonotope'))

% predefine the nonzero uncertainties
up_spec2.p_nonzero = [1 1 zeros(1,dim_u-2)]';
pred = initializeUncertainties(pred,up_spec2);
assert(all(pred.U_init.G == G_u_def(:,1:2),'all'))
assert(isa(pred.U_init,'zonotope'))

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
        assert(isa(pred.U_init,'zonotope'))

        % check number of identified uncertainties
        if strcmp(up_method,"outRandom*")
            assert(size(pred.U_init.G,2) == 2*numUnc)
        else
            assert(size(pred.U_init.G,2) == numUnc)
        end

        if numUnc == dim_u && ~strcmp(up_method,"outRandom*")
            % check if same results for identifying all uncertainties
            assert(all(withinTol(pred_all.U_init.Z,pred.U_init.Z,1e-6),'all'))
        end
    end
end

res = true;
end

% ------------------------------ END OF CODE ------------------------------
