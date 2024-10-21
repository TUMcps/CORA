function [fitness,gp,theta,ypredtrain,fitnessTest,ypredtest,pvals,r2train,r2test,r2val,geneOutputs,geneOutputsTest,geneOutputsVal]=regressmulti_fitfun(evalstr,gp)
%REGRESSMULTI_FITFUN Fitness function for multigene symbolic regression.
%
%   This is the default fitness function for multigene symbolic regression
%   in GPTIPS.
%
%   [FITNESS,GP] = REGRESSMULTI_FITFUN(EVALSTR,GP) returns the FITNESS of
%   the symbolic expression(s) in the cell array EVALSTR using information
%   contained in the GP data struct. Here, FITNESS is the root mean squared
%   prediction error (RMSE) on the training data set.
%
%   [FITNESS,GP,THETA,YPREDTRAIN,FITNESS_TEST,YPREDTEST,PVALS,R2TRAIN,R2TEST,R2VAL]
%   = REGRESSMULTI_FITFUN(EVALSTR,GP) may be used post-run to return the
%   gene coefficients THETA, the prediction of the model on the training
%   data YPREDTRAIN, the RMSE fitness value FITNESS_TEST on the test data
%   set, the prediction of the model on the test data YPREDTEST, the
%   statistical p-values for bias and model terms are returned as PVALS
%   (PVALS only computed if the Statistics Toolbox is present, otherwise an
%   empty variable is returned). Additionally, coefficients of
%   determination (R^2) are returned as R2TRAIN, R2TEST and R2VAL.
%
%   Remarks:
%
%   Each observation of the response variable y is assumed to be an unknown
%   non-linear function of the corresponding observations of the predictor
%   variables x1,..xn.
%
%   Training data:
%
%   The GPTIPS configuration file should populate the following required
%   fields for the training data assuming 'Ntrain' observations on the
%   input and output data. GP.USERDATA.XTRAIN should be a (Ntrain X n)
%   matrix where the ith column contains the Ntrain observations of the ith
%   input variable xi. GP.USERDATA.YTRAIN should be a (Ntrain x 1) vector
%   containing the corresponding observations of the response variable y.
%
%   Testing data:
%
%   The following fields are optional and may be used, post-run, to see how
%   well evolved models generalise to an unseen test data set with Ntest
%   observations. They do not affect the model building process.
%   GP.USERDATA.XTEST should be a (Ntest X n) matrix where the ith column
%   contains the Ntest observations of the ith input variable xi.
%   GP.USERDATA.YTEST should be a (Ntest x 1) vector containing the
%   corresponding observations of the response variable y.
%
%   How multigene symbolic regression works:
%
%   In multigene symbolic regression, each prediction of y is formed by the
%   weighted output of each of the trees/genes in the multigene individual
%   plus a bias term. The number (M) and structure of the trees is evolved
%   automatically during a GPTIPS run (subject to user defined
%   constraints).
%
%   i.e. ypredtrain = c0 + c1*tree1 + ... + cM*treeM
%
%   where c0 = bias term
%         c1,..,cM are the weights
%         M is the number of genes/trees comprising the current individual
%
%   The weights (i.e. regression coefficients) are automatically determined
%   by a least squares procedure for each multigene individual and are
%   stored in GP.FITNESS.RETURNVALUES for future use.
%
%   Remarks:
%
%   Because the GP structure is modified within this function (i.e. the
%   field GP.FITNESS.RETURNVALUES is used to store the computed weighting
%   coefficients for each gene) the GP structure must be returned as an
%   output argument.
%
%   This fitness function is used for multigene symbolic regression for
%   GPDEMO2, GPDEMO3 and GPDEMO4 (the configuration files for these are
%   GPDEMO2_CONFIG.M and GPDEMO3_CONFIG.M respectively) but it can and
%   should be used for the user's own non-linear regression problems.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also REGRESSMULTI_FITFUN_VALIDATE, GPDEMO2_CONFIG, GPDEMO3_CONFIG,
%   GPDEMO4_CONFIG, GPDEMO2, GPDEMO3

%defaults in case of early exit
theta=[];ypredtrain=[];fitnessTest=[];ypredtest=[];pvals=[];
r2train=[];r2test=[];r2val=[];geneOutputs=[];geneOutputsTest=[];
geneOutputsVal=[];

% process evalstr with regex to allow direct access to data matrices
pat = 'x(\d+)';
evalstr = regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)');
y = gp.userdata.ytrain;
numData = gp.userdata.numytrain;
numGenes = numel(evalstr);

%set up a matrix to store the tree outputs plus a bias column of ones
geneOutputs = ones(numData,numGenes+1);

%eval each gene in the current individual
for i = 1:numGenes
    ind = i + 1;
    eval(['geneOutputs(:,ind)=' evalstr{i} ';']);
    
    %check for nonsensical answers and break out early with an 'inf' if so
    if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
        fitness = Inf;
        gp.fitness.returnvalues{gp.state.current_individual} = [];
        return
    end
end

%only calc. weighting coeffs during an actual run or if forced
if ~gp.state.run_completed || gp.state.force_compute_theta
    
    %set gp.userdata.bootSample to true to resample data for weights computation
    
    %prepare LS matrix
    if gp.userdata.bootSample
        sampleInds = bootsample(geneOutputs,gp.userdata.bootSampleSize);
        goptrans = geneOutputs(sampleInds,:)';
        prj = goptrans * geneOutputs(sampleInds,:);
        ysample = y(sampleInds);
    else
        goptrans = geneOutputs';
        prj = goptrans * geneOutputs;
    end
    
    %calculate tree weight coeffs using SVD based least squares
    %normal equation
    try
        if gp.userdata.bootSample
            theta = pinv(prj) * goptrans * ysample;
        else
            theta = pinv(prj) * goptrans * y;
        end
    catch
        theta = [];
        fitness = Inf;
        gp.fitness.returnvalues{gp.state.current_individual} = [];
        return;
    end
    
    %assign bad fitness if any coeffs NaN or Inf
    if any(isinf(theta)) || any(isnan(theta))
        theta = [];
        fitness = Inf;
        gp.fitness.returnvalues{gp.state.current_individual} = [];
        return;
    end
    
    %write coeffs to returnvalues field for storage
    gp.fitness.returnvalues{gp.state.current_individual} = theta;
    
else %if post-run, get stored coeffs from return value field
    theta = gp.fitness.returnvalues{gp.state.current_individual};
end

%calc. prediction of full training data set using the estimated weights
ypredtrain = geneOutputs * theta;

%calculate RMS prediction error (fitness)
err = gp.userdata.ytrain - ypredtrain;
fitness = sqrt(((err'*err)/numData));

%--below is for post-run evaluation of models, it is not used during a GPTIPS run--

if gp.state.run_completed
    
    %compute r2 for training data
    r2train = 1 - sum( (gp.userdata.ytrain-ypredtrain).^2 )/sum( (gp.userdata.ytrain-mean(gp.userdata.ytrain)).^2 );
    plotValidation = 0;
    
    %process validation data if present
    if (isfield(gp.userdata,'xval')) && (isfield(gp.userdata,'yval')) && ...
            ~isempty(gp.userdata.xval) && ~isempty(gp.userdata.yval)
        
        plotValidation = 1;
        evalstr = strrep(evalstr,'.xtrain','.xval');
        yval = gp.userdata.yval;
        numData = length(yval);
        
        %set up a matrix to store the tree outputs plus a bias column of ones
        geneOutputsVal = zeros(numData,numGenes + 1);
        geneOutputsVal(:,1) = ones;
        
        %eval each tree
        for i=1:numGenes
            ind = i+1;
            eval(['geneOutputsVal(:,ind)=' evalstr{i} ';']);
        end
        
        ypredval = geneOutputsVal*theta; %create the prediction  on the validation data
        fitness_val = sqrt(mean((gp.userdata.yval - ypredval).^2));
        
        %compute r2 for validation data
        r2val = 1 - sum( (gp.userdata.yval - ypredval).^2 )/sum( (gp.userdata.yval - mean(gp.userdata.yval)).^2 );
        evalstr = strrep(evalstr,'.xval','.xtrain');
    else
        r2val = [];
    end %end of validation data calcs
    
    
    %process test data if present
    plotTest = 0;
    if (isfield(gp.userdata,'xtest')) && (isfield(gp.userdata,'ytest')) && ...
            ~isempty(gp.userdata.xtest) && ~isempty(gp.userdata.ytest)
        
        plotTest = 1;
        evalstr = strrep(evalstr,'.xtrain','.xtest');
        ytest = gp.userdata.ytest;
        numData = length(ytest);
        
        %set up a matrix to store the tree outputs plus a bias column of ones
        geneOutputsTest = zeros(numData,numGenes+1);
        geneOutputsTest(:,1) = ones;
        
        %eval each tree
        for i=1:numGenes
            ind = i + 1;
            eval(['geneOutputsTest(:,ind)=' evalstr{i} ';']);
        end
        
        ypredtest = geneOutputsTest * theta; %create the prediction on the testing data
        fitnessTest = sqrt(mean((gp.userdata.ytest - ypredtest).^2));
        
        %compute r2 for test data
        r2test = 1 - sum( (gp.userdata.ytest - ypredtest).^2 )/sum( (gp.userdata.ytest - mean(gp.userdata.ytest)).^2 );
    end
    
    
    %calc statistical analysis of gene significance on training data
    %(if stats toolbox is present)
    if gp.userdata.stats && gp.info.toolbox.stats
        % Regress tree outputs (and bias) against y train data and get stats
        wstate = warning;warning off;
        i_max = min([size(y,1),10000]);
        stats = regstats(y(1:i_max,:),geneOutputs(1:i_max,2:end));
        warning(wstate);
        pvals = stats.tstat.pval;
    else
        pvals = [];
    end
end

%if graphs required
if gp.state.run_completed && gp.userdata.showgraphs
    
    if ~isempty(gp.userdata.name)
        setname = ['Data: ' gp.userdata.name];
    else
        setname='';
    end
    i_max = min([size(y,1),500]);
    
    %model predictions plots
    figure('name','GPTIPS 2 Multigene regression. Model predictions.','numbertitle','off');
    subplot(1+plotTest+plotValidation,1,1);
    plot(ypredtrain(1:i_max,:),'Color',[0.85 0.33 0.1]);
    hold on;
    plot(gp.userdata.ytrain(1:i_max,:),'Color',[0 0.45 0.74]);
    axis tight;
    ylabel('y');
    xlabel('Data point');
    legend('Predicted','Actual');
    title({setname,...
        ['RMS training set error: ' num2str(fitness) '  R^2: ' num2str(r2train)]},'interpreter','tex');
    hold off;
    
    if plotTest
        subplot(2+plotValidation,1,2);
        plot(ypredtest(1:i_max,:),'Color',[0.85 0.33 0.1]);
        hold on;
        plot(gp.userdata.ytest(1:i_max,:),'Color',[0 0.45 0.74]);
        axis tight;
        ylabel('y');
        xlabel('Data point');
        title(['RMS test set error: ' num2str(fitnessTest) '  R^2: ' num2str(r2test)],'interpreter','tex');
        hold off
    end
    
    if plotValidation
        subplot(2+plotValidation,1,3);
        plot(ypredval(1:i_max,:),'Color',[0.85 0.33 0.1]);
        hold on;
        plot(gp.userdata.yval(1:i_max,:),'Color',[0 0.45 0.74]);
        axis tight;
        ylabel('y');
        xlabel('Data point');
        title(['RMS validation set error: ' num2str(fitness_val) '  R^2: ' num2str(r2val)],'interpreter','tex');
        hold off
    end
    
    %scatterplots
    scatterFig = figure('name','GPTIPS 2 Multigene regression. Model prediction scatterplot.','numbertitle','off');
    subplot(1+plotTest+plotValidation,1,1);
    minval = min([gp.userdata.ytrain;ypredtrain]);
    maxval = max([gp.userdata.ytrain;ypredtrain]);
    axis([minval maxval minval maxval]);
    ilineTr = line([minval maxval], [minval maxval]); 
    set(ilineTr,'color','black','LineWidth',1);hold on;
    scatter(gp.userdata.ytrain,ypredtrain,'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','none');
    box on;grid on;hold off;
    ylabel('Predicted');
    xlabel('Actual');
    title({setname,['RMS training set error: ' num2str(fitness) '  R^2: ' num2str(r2train)]},'interpreter','tex');
    
    %add scatter plot for test data, if present
    if plotTest
        subplot(2+plotValidation,1,2);
        minval = min([gp.userdata.ytest;ypredtest]);
        maxval = max([gp.userdata.ytest;ypredtest]);
        axis([minval maxval minval maxval]);
        ilineTest = line([minval maxval], [minval maxval]);
        set(ilineTest,'color','black','LineWidth',1);hold on;
        scatter(gp.userdata.ytest,ypredtest,'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','none');
        box on;grid on;hold off;
        ylabel('Predicted');
        xlabel('Actual');
        title(['RMS test set error: ' num2str(fitnessTest) '  R^2: ' num2str(r2test)],'interpreter','tex');
    end
    
    %add scatter plot for validation data, if present
    if plotValidation
        figure(scatterFig);
        subplot(2+plotTest,1,2+plotTest);
        minval = min([gp.userdata.yval;ypredval]);
        maxval = max([gp.userdata.yval;ypredval]);
        axis([minval maxval minval maxval]);
        ilineVal = line([minval maxval], [minval maxval]);
        set(ilineVal,'color','black','LineWidth',1);hold on;
        scatter(gp.userdata.yval,ypredval,'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','none');
        box on;grid on;hold off;
        title(['RMS validation set error: ' num2str(fitness_val) '  R^2: ' num2str(r2val)],'interpreter','tex');
        ylabel('Predicted');xlabel('Actual');
    end
    
    %gene weights & significance
    if  gp.info.toolbox.stats && gp.userdata.stats
        
        %generate x labels for bar graphs
        geneLabels = {'Bias'};
        for i = 1:numGenes
            geneLabels{i+1} = ['Gene ' int2str(i)];
        end
        
        %plot gene weights and offset
        statFig = figure; coeffsAx = subplot(2,1,1);
        set(statFig,'name','GPTIPS 2 P-values of model genes (on training data)','numbertitle','off');
        geneBar = bar(coeffsAx,stats.beta);
        set(coeffsAx,'xtick',1:(numGenes+1));
        set(coeffsAx,'xticklabel',geneLabels);
        title(coeffsAx,{setname,'Gene weights'});
        
        %plot p-vals
        pvalsAx = subplot(2,1,2);
        pvalBar = bar(pvalsAx,stats.tstat.pval);
        
        if ~verLessThan('matlab','8.4') %R2014b
            pvalBar.FaceColor = [0 0.45 0.74];
            pvalBar.BaseLine.Visible = 'off';
            geneBar.FaceColor = [0 0.45 0.74];
        else
            pvalsBars = get(pvalsAx,'Children');
            set(pvalsBars,'FaceColor',[0 0.45 0.74],'ShowBaseLine','off');
            coeffsBars = get(coeffsAx,'Children');
            set(coeffsBars,'FaceColor',[0 0.45 0.74]);
        end
        
        set(pvalsAx,'xtick',1:(numGenes+1));
        set(pvalsAx,'xticklabel',geneLabels);
        title(pvalsAx,'P value');xlabel(['R^2 = ' num2str(r2train)]);
        
    end
end