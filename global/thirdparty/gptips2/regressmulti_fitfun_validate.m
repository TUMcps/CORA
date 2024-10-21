function [valfitness,gp,ypredval] = regressmulti_fitfun_validate(gp)
%REGRESSMULTI_FITFUN_VALIDATE Evaluate current 'best' multigene regression model on validation data set.
%
%   [FITNESS,GP,YPREDVAL] = REGRESSMULTI_FITFUN_VALIDATE(GP) returns the
%   FITNESS of the the current 'best' multigene model contained in the GP
%   data structure (as evaluated on the training data). FITNESS is the root
%   mean squared prediction (RMSE) error on the validation data set.
%
%   Remarks on the use of validation data:
%
%   This is used in conjunction with the multigene regression fitness
%   function REGRESSMULTI_FITFUN.
%
%   This function is called at the end of each generation by adding the
%   following line to the user's run configuration file:
%
%   GP.USERDATA.USER_FCN = @REGRESSMULTI_FITFUN_VALIDATE
%
%   In addition, the user's GPTIPS configuration file should populate the
%   following required fields for the training data assuming Nval
%   observations on the input and output data: GP.USERDATA.XVAL should be a
%   (Nval x n) matrix where the ith column contains the Nval observations
%   of the ith input variable xi. GP.USERDATA.YVAL should be a (Nval x 1)
%   vector containing the corresponding observations of the response
%   variable y.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also REGRESSMULTI_FITFUN, GP_USERFCN

% on first call, check that validation data is present
if (gp.state.count == 1) && (~isfield(gp.userdata,'xval')) || (~isfield(gp.userdata,'yval')) || ...
        isempty(gp.userdata.xval) || isempty(gp.userdata.yval)
    error('Cannot perform holdout validation because no validation data was found.');
end

%get evalstr for current 'best' individual on training data
evalstr = gp.results.best.eval_individual;

%process evalstr with regex to allow direct access to data
pat = 'x(\d+)';
evalstr = regexprep(evalstr,pat,'gp.userdata.xval(:,$1)');
y = gp.userdata.yval;

numDataPoints = size(y,1);
numGenes = length(evalstr);

%set up a matrix to store the tree outputs plus a bias column of ones
gene_outputs = ones(numDataPoints,numGenes+1);

%eval each gene in the current individual
for i=1:numGenes
    ind = i + 1;
    eval(['gene_outputs(:,ind)=' evalstr{i} ';']);
    
    %check for nonsensical answers and break out early with an 'inf' if so
    if  any(~isfinite(gene_outputs(:,ind))) || any(~isreal(gene_outputs(:,ind)))
        valfitness = Inf;
        return
    end
end

%retrieve regression weights
theta = gp.results.best.returnvalues;

%calc. prediction of validation data set using the retreived weights
try
    ypredval = gene_outputs * theta;
    %calculate RMS prediction error on validation data
    valfitness = sqrt(mean((y-ypredval).^2));
catch ME
    display(ME)
    valfitness = Inf;
end


%on 1st gen, initialise validation set info in the GP structure
if gp.state.count == 1
    gp.results.history.valfitness(1:gp.runcontrol.num_gen,1) = 0;
    gp.results.history.valfitness(1) = valfitness;
    gp.results.valbest = gp.results.best;
    gp.results.valbest.valfitness = valfitness;
end

gp.results.best.valfitness = valfitness;
gp.results.history.valfitness(gp.state.count,1) = valfitness;

%update 'best' validation if fitness is better or its the same with less
%nodes/complexity
if gp.fitness.complexityMeasure
    bcomp = gp.results.best.complexity;
    vcomp = gp.results.valbest.complexity;
else
    bcomp = gp.results.best.nodecount;
    vcomp = gp.results.valbest.nodecount;
end

if valfitness < gp.results.valbest.valfitness || ...
        (valfitness == gp.results.valbest.valfitness && bcomp < vcomp)
    gp.results.valbest = gp.results.best;
    gp.results.valbest.valfitness = valfitness;
    gp.results.valbest.foundatgen = gp.state.count - 1;
end