function gpmodel = gpmodel2struct(gp,ID,tbxStats,createSyms,modelStruc,fastSymMode)
%GPMODEL2STRUCT Create a struct describing a multigene regression model.
%
%   GPMODEL = GPMODEL2STRUCT(GP,ID) gets the multigene regression model
%   specified by the identifier ID from the GPTIPS struct GP and returns
%   model info and performance data as a struct GPMODEL.
%
%   The GPMODEL2STRUCT struct contains a variety of model information
%   including symbolic  object versions of the multigene regression  model
%   and its constituent genes as well as model performance metrics for the
%   training, validation and test data sets (if present). This struct - as
%   well as being a convenient store of model performance data - may also
%   be used as an input to several other GPTIPS model analysis functions,
%   e.g.
%
%   GPPRETTY(GP,GPMODEL)
%
%   GPMODELREPORT(GP,GPMODEL)
%
%   RUNTREE(GP,GPMODEL)
%
%   GPMODEL2SYM(GP,GPMODEL)
%
%   GPMODEL2MFILE(GP,GPMODEL)
%
%   GPMODEL2FUNC(GP,GPMODEL)
%
%   GPMODELGENES2MFILE(GP,GPMODEL)
%
%   DRAWTREES(GP,GPMODEL)
%
%   Note:
%
%   The GPMODEL struct is functionally identical to that returned by the
%   function GPGENES2MODEL wherein a 'new' multigene regression model may
%   be constructed from the unique genes in a population.
%
%   Additionally:
%
%   GPMODEL = GPMODEL2STRUCT(GP,'best') gets info on the 'best' model in
%   the population (as evaluated on training data).
%
%   GPMODEL = GPMODEL2STRUCT(GP,'valbest') gets info on the 'best' model
%   (as evaluated on the validation data set (if this data exists).
%
%   GPMODEL = GPMODEL2STRUCT(GP,'testbest') gets info on the 'best' model
%   (as evaluated on the validation data set (if this data exists).
%
%   GPMODEL = GPMODEL2STRUCT(GP,ID,TBXSTATS) where TBXSTATS is TRUE does
%   the same but computes additional model performance stats (on the
%   training data) using the Statistics Toolbox (TBXSTATS default is
%   FALSE). These may be found in the field GPMODEL.TRAIN.TBXSTATS which
%   contains the output of the Statistics Toolbox REGSTATS function.
%
%   GPMODEL = GPMODEL2STRUCT(GP,ID,TBXSTATS,CREATESYMS) where CREATESYMS is
%   FALSE does the same but does not generate any symbolic math objects
%   using the Symbolic Math Toolbox. (CREATESYMS default is TRUE).
%
%   Remarks:
%
%   This function always returns a GPMODEL struct even if the supplied
%   model identifier is invalid (e.g. if 'valbest' is specified but there
%   is no validation data) or if the model is invalid in some other way
%   (e.g. there were no gene weights computed for it due to non-finite gene
%   outputs on the training data). The field GPMODEL.VALID in the returned
%   struct is TRUE for valid models and FALSE otherwise. If the model is
%   invalid then the reason for the model's invalidity can always be found
%   in the field GPMODEL.INVALIDREASON. E.g. 'Invalid model index
%   supplied.'
%
%   Furthermore - for a valid model - if there was 'test' set data supplied
%   during the GP run then the field GPMODEL.TEST.WARNING will be TRUE if
%   there was a problem predicting the test set values amd the reason for
%   the problem will be in the field GPMODEL.TEST.WARNINGREASON (e.g.
%   'Non-finite or non-real predictions on testing data.'). Similarly for
%   validation data, if there was a problem predicting it then the field
%   GPMODEL.VAL.WARNING is TRUE and the reason found in
%   GPMODEL.VAL.WARNINGREASON.
%
%   It is quite possible for a model to be formally 'valid' (i.e. predicts
%   OK on the training data) but fails on the test or validation data (e.g.
%   due to division by zero etc.)
%
%   This function computes additional model performance stats using the
%   sub-function REGRESSMULTI_FITFUN_FULL_STATS. So, if you change the
%   multigene regression fitness function REGRESSMULTI_FITFUN then you will
%   also need to change REGRESSMULTI_FITFUN_FULL_STATS.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPMODEL2MFILE, GPMODEL2FUNC, GENES2GPMODEL, GPMODEL2SYM,
%            GPMODELREPORT, DRAWTREES, RUNTREE, GPMODELGENES2MFILE

if nargin < 2
    gpmodel.valid = false;
    gpmodel.invalidReason = 'The function GPMODEL2STRUCT requires at least 2 arguments, e.g. GPMODEL2STRUCT(GP,''BEST'')';
    return;
end

if ~strncmpi('regressmulti',func2str(gp.fitness.fitfun),12)
    error('GPMODEL2STRUCT may only be used to extract model data from a GP structure with a population containing multigene symbolic regression models.');
end

if nargin < 6 || isempty(fastSymMode)
    fastSymMode = false;
end

%scan model genes for input frequency, depth, complexity? (default = yes).
%(setting to false is a bit quicker but doesn't compute the structural info)
if nargin < 5 || isempty(modelStruc)
    modelStruc = true;
end

%create symbolic objects of overall model and genes? (default = yes).
%(setting to false is significantly quicker but you don't get the symbolic
%math object for each gene)
if nargin < 4 || isempty(createSyms)
    createSyms = true;
end

%compute statistics toolbox stats? (default = no).
if nargin < 3 || isempty(tbxStats)
    tbxStats = false;
end

%set up default test, train and validation values
gpmodel.train.r2 = 0;
gpmodel.train.rmse = Inf;
gpmodel.train.mse = Inf;
gpmodel.train.sse = Inf;
gpmodel.train.mae = Inf;
gpmodel.train.maxe = Inf;
gpmodel.train.err = [];
gpmodel.train.ypred = [];
gpmodel.train.gene_outputs = [];
gpmodel.train.datapoints = [];
gpmodel.train.warning = false;
gpmodel.train.warningReason = '';

gpmodel.val.r2 = 0;
gpmodel.val.rmse = Inf;
gpmodel.val.mse = Inf;
gpmodel.val.sse = Inf;
gpmodel.val.mae = Inf;
gpmodel.val.maxe = Inf;
gpmodel.val.err = [];
gpmodel.val.ypred = [];
gpmodel.val.gene_outputs = [];
gpmodel.val.datapoints = [];
gpmodel.val.warning = false;
gpmodel.val.warningReason = '';

gpmodel.test.r2 = 0;
gpmodel.test.rmse = Inf;
gpmodel.test.mse = Inf;
gpmodel.test.sse = Inf;
gpmodel.test.mae = Inf;
gpmodel.test.maxe = Inf;
gpmodel.test.err = [];
gpmodel.test.ypred = [];
gpmodel.test.gene_outputs = [];
gpmodel.test.datapoints = [];
gpmodel.test.warning = false;
gpmodel.test.warningReason = '';

%encoded genes supplied in cell array
cellgenes = false;

%parse user input for supplied numerical model identifier
if isnumeric(ID) && numel(ID) == 1
    
    if  ~mod(ID,1) && ID > 0 && ID <= gp.runcontrol.pop_size
        
        %get encoded trees, eval trees and return values
        treestrs = gp.pop{ID};
        evaltreestrs = tree2evalstr(gp.pop{ID},gp);
        rtnVals = gp.fitness.returnvalues{ID};
    else
        gpmodel.valid = false;
        gpmodel.invalidReason = 'Invalid model index supplied.';
        return;
    end
    
elseif ischar(ID) && strcmpi(ID,'best') %or 'best' on training
    
    treestrs = gp.results.best.individual;
    evaltreestrs = gp.results.best.eval_individual;
    rtnVals = gp.results.best.returnvalues;
    
elseif  ischar(ID) && strcmpi(ID,'valbest') %or 'best' on validation data
    
    % check that validation data is present
    if ~isfield(gp.results,'valbest')
        gpmodel.valid = false;
        gpmodel.invalidReason ='No validation data was found.';
        return;
    end
    
    treestrs = gp.results.valbest.individual;
    evaltreestrs = gp.results.valbest.eval_individual;
    rtnVals = gp.results.valbest.returnvalues;
    
elseif  ischar(ID) && strcmpi(ID,'testbest') %or 'best' on test data
    
    % check that test data is present
    if ~isfield(gp.results,'testbest')
        gpmodel.valid = false;
        gpmodel.invalidReason ='No test data was found.';
        return;
    end
    
    treestrs = gp.results.testbest.individual;
    evaltreestrs = gp.results.testbest.eval_individual;
    rtnVals = gp.results.testbest.returnvalues;
    
    
    %if user supplied list of genes in encoded form,e.g. {'n(x1,f(x8),x11)','c(x2,x8)'}
elseif iscell(ID)
    
    treestrs = ID;
    evaltreestrs = tree2evalstr(treestrs,gp);
    cellgenes = true;
    
    %return values need to be computed by the fitness function
    gp.state.force_compute_theta = true;
    gp.state.run_completed = false;
    gp.userdata.showgraphs = false;
    [fitness,gp,rtnVals] = feval(gp.fitness.fitfun,evaltreestrs,gp);
    
    if isinf(fitness)
        gpmodel.valid = false;
        gpmodel.invalidReason ='Supplied gene list gave non-finite values for predictions of training data';
        return
    end
    
elseif isstruct(ID)
    
    gpmodel.valid = false;
    gpmodel.invalidReason ='Models already in struct format unsupported by gpmodel2struct';
    return
    
else %otherwise user did not supply a valid multigene model selector
    
    gpmodel.valid = false;
    gpmodel.invalidReason ='Invalid model selector supplied.';
    return
end

%get multigene regression stats
wstate = warning; warning off;
gpmodel = regressmulti_fitfun_full_stats(gpmodel,evaltreestrs,gp,rtnVals,tbxStats);
warning(wstate);

%compute some tree structural information
if modelStruc
    
    gpmodel.numNodes = getnumnodes(treestrs);
    gpmodel.expComplexity = getcomplexity(treestrs);
    inputVec = gpmodelvars(gp,ID);
    
    %inputs
    count = 0;
    for i=1:numel(inputVec)
        if inputVec(i) > 0
            count = count + 1;
            if ~isempty(gp.nodes.inputs) && i <= numel(gp.nodes.inputs.names) && ~isempty(gp.nodes.inputs.names{i})
                gpmodel.inputs{count} = strtrim(gp.nodes.inputs.names{i});
            else
                gpmodel.inputs{count} = ['x' num2str(i)];
            end
        end
    end
    gpmodel.numInputs = numel(gpmodel.inputs);
    gpmodel.output = gp.nodes.output.name;
    
    %max depth
    maxDepth = 1;
    for i=1:numel(treestrs)
        maxDepth = max(maxDepth,getdepth(treestrs{i}));
    end
    gpmodel.maxDepth = maxDepth;
    
end

%create sym objects
if createSyms && gp.info.toolbox.symbolic
    
    if cellgenes
        ID = horzcat(ID,rtnVals);
    end
    
    [gpmodel.sym,gpmodel.genes.geneSyms] = gpmodel2sym(gp,ID,fastSymMode,true);
else
    gpmodel.sym = [];
    gpmodel.genes.geneSyms = [];
end

%original encoded genes
gpmodel.about = 'A struct representing a multigene regression model.';
gpmodel.genes.geneStrs = treestrs;
gpmodel.genes.geneWeights = rtnVals;
gpmodel = orderfields(gpmodel);
gpmodel.source = 'gpmodel2struct';

function gpmodel = regressmulti_fitfun_full_stats(gpmodel, evalstr,gp,theta,tbxStats)
%REGRESSMULTI_FITFUN_FULL_STATS updates GPMODEL struct with stats for an existing multigene regression model.
%
%    GPMODEL = REGRESSMULTI_FITFUN_FULL_STATS(EVALSTR,GP,THETA,GRAPHS,TBXST
%    ATS) updates a structure GPMODEL containing performance data about the
%    multigene symbolic regression individual represented by EVALSTR and
%    THETA.
%
%    Remarks:
%
%    Utility function. Reproduces functionality of REGRESSMULTI_FITFUN in a
%    more useful (but slower) form for offline use. It calculates a variety
%    of performance statistics which can be accessed via the fields of
%    GPMODEL. Also sets warning flags for 'test' and 'validation' data sets
%    if the model fails to predict for them.

%statistics toolbox stats off by default
if nargin < 4
    tbxStats = false;
end

%ensure gene weights in matrix form
if iscell(theta)
    theta = theta{1};
end

gpmodel.valid = true;
gpmodel.invalidReason = '';

%if the gene weights vector is empty it means that during the
%fitness evaluation it could not be generated due to either
%non-finite or complex model outputs OR numerical problem in the
%SVD least squares computation of weights. This is a model "showstopper"
%and so the model is marked as 'invalid'
if numel(evalstr) ~= (numel(theta)-1)
    gpmodel.valid = false;
    gpmodel.invalidReason = 'Empty gene weights vector: probably because model output contained complex or non-finite values on training data.';
    return
end

%process evalstr with regex to allow direct access to data matrices
pat = 'x(\d+)';
evalstr = regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)');
numTrainData = numel(gp.userdata.ytrain);
numGenes = length(evalstr);
gpmodel.genes.num_genes = numGenes;

%set up a matrix to store the tree outputs plus a bias column of ones
geneOutputsTrain = ones(numTrainData,numGenes+1);

%eval each gene in the current individual on training data
for i=1:numGenes
    ind = i + 1;
    eval(['geneOutputsTrain(:,ind)=' evalstr{i} ';']);
end

%add raw gene outputs to struct (not bias term)
gpmodel.train.gene_outputs = geneOutputsTrain(:,2:end);
gpmodel.train.datapoints = numTrainData;

%check for nonsensical answers - exit with invalid model if any found
if  any(any(~isfinite(geneOutputsTrain))) || any(any(~isreal(geneOutputsTrain)))
    gpmodel.valid = false;
    gpmodel.invalidReason = 'Non-finite or non-real gene output values on training data.';
    return
end

%calc contribution of weighted individual genes on training data
gpmodel.train.ypred_genes = geneOutputsTrain.*repmat(theta',numTrainData,1);

%calc overall prediction of training data
gpmodel.train.ypred = geneOutputsTrain*theta;

%error (training data)
errTrain = gp.userdata.ytrain - gpmodel.train.ypred;
gpmodel.train.err = errTrain;

%SSE (training data)
gpmodel.train.sse = (errTrain'*errTrain);

%MSE (training data)
gpmodel.train.mse = gpmodel.train.sse/numTrainData;

%RMS prediction error - training data
gpmodel.train.rmse = sqrt(gpmodel.train.mse);

%r2 for training data
r2train = 1 - (gpmodel.train.sse /sum( (gp.userdata.ytrain-mean(gp.userdata.ytrain)).^2 ) );
gpmodel.train.r2 = r2train;

%mean absolute error (training)
gpmodel.train.mae = mean(abs(errTrain));

%max abs error (training)
gpmodel.train.maxe = max(abs(errTrain));

%process validation data if present
if isfield(gp.results,'valbest')
    
    evalstr = strrep(evalstr,'.xtrain','.xval');
    numValData = length(gp.userdata.yval);
    
    geneOutputsVal = zeros(numValData,numGenes+1);
    geneOutputsVal(:,1) = ones;
    
    for i=1:numGenes
        ind = i + 1;
        eval(['geneOutputsVal(:,ind)=' evalstr{i} ';']);
    end
    
    gpmodel.val.gene_outputs = geneOutputsVal(:,2:end);
    gpmodel.val.datapoints = numValData;
    
    %flag warning if non-real or complex predictions on validation data
    if  any(any(~isfinite(geneOutputsVal))) || any(any(~isreal(geneOutputsVal)))
        gpmodel.val.rmse = Inf;
        gpmodel.val.warning = true;
        gpmodel.val.warningReason = 'Non-finite or non-real predictions on validation data.';
    end
    
    %compute model stats and predictions if no warning
    if ~gpmodel.val.warning
        
        gpmodel.val.ypred = geneOutputsVal*theta; %create the prediction on the validation data
        
        %error validation data
        errVal = gp.userdata.yval - gpmodel.val.ypred;
        gpmodel.val.err = errVal;
        
        %sse for validation data
        gpmodel.val.sse = errVal'*errVal;
        
        %MSE validation data
        gpmodel.val.mse = gpmodel.val.sse/numValData;
        
        %rmse validation data
        gpmodel.val.rmse = sqrt(gpmodel.val.mse);
        
        %R2 for validation data
        gpmodel.val.r2 = 1 - (gpmodel.val.sse/sum( (gp.userdata.yval - mean(gp.userdata.yval)).^2 ));
        
        %mean absolute error for validation data
        gpmodel.val.mae = mean(abs(errVal));
        
        %max abs error (validation)
        gpmodel.val.maxe = max(abs(errVal));
        
    end
    
    evalstr = strrep(evalstr,'.xval','.xtrain');
    
else
    gpmodel.val.rmse = Inf;
    gpmodel.val.warning = true;
    gpmodel.val.warningReason = 'No validation data was found.';
end %end of validation data calcs


%process test data
if (isfield(gp.userdata,'xtest')) && (isfield(gp.userdata,'ytest')) && ...
        ~isempty(gp.userdata.xtest) && ~isempty(gp.userdata.ytest)
    
    evalstr = strrep(evalstr,'.xtrain','.xtest');
    numTestData = length(gp.userdata.ytest);
    
    geneOutputsTest = zeros(numTestData,numGenes+1);
    geneOutputsTest(:,1) = ones;
    
    for i=1:numGenes
        ind = i + 1;
        eval(['geneOutputsTest(:,ind)=' evalstr{i} ';']);
    end
    
    gpmodel.test.gene_outputs = geneOutputsTest(:,2:end);
    gpmodel.test.datapoints = numTestData;
    
    if  any(any(~isfinite(geneOutputsTest))) || any(any(~isreal(geneOutputsTest)))
        gpmodel.test.rmse = Inf;
        gpmodel.test.warning = true;
        gpmodel.test.warningReason = 'Non-finite or non-real predictions on testing data.';
    end
    
    if ~gpmodel.test.warning
        
        gpmodel.test.ypred = geneOutputsTest*theta;
        
        %error test data
        errTest = gp.userdata.ytest - gpmodel.test.ypred;
        gpmodel.test.err = errTest;
        
        %sse for test data
        gpmodel.test.sse = errTest'*errTest;
        
        %MSE test
        gpmodel.test.mse = gpmodel.test.sse/numTestData;
        
        %RMSE test data
        gpmodel.test.rmse = sqrt(gpmodel.test.mse);
        
        %r2 test data
        gpmodel.test.r2 = 1- (gpmodel.test.sse/sum( (gp.userdata.ytest - mean(gp.userdata.ytest)).^2 ));
        
        %mean absolute error test data
        gpmodel.test.mae = mean(abs(errTest));
        
        %max abs error (testing)
        gpmodel.test.maxe = max(abs(errTest));
        
    end
else
    gpmodel.test.rmse = Inf;
    gpmodel.test.warning = true;
    gpmodel.test.warningReason = 'No test data was found.';   
end

%calc statistical analysis of gene significance & other stats on training data
%(if stats toolbox is present)
if tbxStats && gp.info.toolbox.stats
    stats = regstats(gp.userdata.ytrain,geneOutputsTrain(:,2:end));
    gpmodel.train.pvals = stats.tstat.pval;
    gpmodel.train.tbxStats = stats;
else
    gpmodel.train.pvals = [];
    gpmodel.train.tbxStats = [];
end


