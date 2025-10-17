function [conservatism, coverage, Y] = evaluateData(pred,data)
% evaluateData - evaluate predictor over data
%
% Syntax:
%    [conservatism, coverage, Y] = evaluateData(pred,data)
%
% Inputs:
%    pred - confPredictor object
%    data - evaluation data, either specified as trajectory objects (where
%           traj(i).x is the i'th input and traj(i).y the corresponding 
%           output) or as struct with data.input(:,i) and data.target(:,i)
%
% Outputs:
%    conservatism - mean prediction set size
%    coverage - ratio of contained data points to number of data points
%    Y - cell array of prediction sets for all data points
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

predictionSetSize_sum = 0;
numsCorrect = 0;

if ~isa(data,'trajectory')
    % transform to trajectory objects
    traj = [];
    u = zeros(pred.sys.nrOfInputs,1);
    for i = 1:size(data.input,2)
        traj = [traj; trajectory(u,data.input(:,i),data.target(:,i),[],1)];
    end
    data = traj;
end

% compute predicted output sets for the calibration data
n_data = length(data);
Y = cell(n_data,1);
if (isa(pred.U,'zonotope') && ~isfinite(sum(pred.U.G,'all'))) || ...
        (isa(pred.U,'interval') && (~isfinite(sum(pred.U.sup,'all')) || ~isfinite(sum(pred.U.inf,'all'))))
    % no valid uncertainty set
    conservatism = inf;
    coverage = 0;
    return
end
for i=1:n_data
    x_i = traj(i).x;
    Y{i} = predict(pred,x_i);
end

% evaluate conservatism and coverage
if strcmp(pred.task, "class")
    % classification: compute number of classes predicted
    vec = 1: length(traj(1).y);
    for i=1:n_data
        class_true = vec(logical(traj(i).y));

        if strcmp(pred.type,'conformalAdd')
            % Y{i} is a vector of predicted classes
            [num, numCorrect] = aux_computeNumClassesCP(Y{i},class_true);
        else
            % Y{i} is a prediction set
            [num, numCorrect] = aux_computeNumClasses(Y{i},class_true);
        end
        predictionSetSize_sum = predictionSetSize_sum + num;
        numsCorrect = numsCorrect + numCorrect;
    end
else
    % regression: compute size of the predicted output set
    for i=1:n_data
        if contains(Y{i},traj(i).y)
            numsCorrect=numsCorrect + 1;
        else
        end
        if dim(Y{i}) > 2
            Yred = reduce(Y{i},'girard',10);
        else
            Yred = Y{i};
        end
        predictionSetSize_sum = predictionSetSize_sum + aux_conservatism(Yred);
    end
end
conservatism = predictionSetSize_sum / n_data;
coverage = numsCorrect / n_data;
end


% Auxiliary functions -----------------------------------------------------


function c = aux_conservatism(R)
if R.dim <= 5
    c = volume(R);
else
    % c = volume(R, 'reduce');
    c = 0;
    for i=0:ceil(R.dim/3)-1
        idx_end = min(3*i+3,dim(R));
        c = c + volume(project(R,3*i+1:idx_end));
    end
end
end

function [num, numCorrect] = aux_computeNumClasses(Y,class_true)
% compute the number of classes which are feasible for a given output set

optionsOptim = optimoptions('linprog','Algorithm', 'interior-point', 'Display','off');
if isa(Y,'interval')
    Y = zonotope(Y);
end
eta = size(generators(Y),2);
n_y = dim(Y);
num = 0;
numCorrect = 0;
for ii = 1:n_y
    % test for each class ii
    if eta == 0 || sum(abs(generators(Y)),'all') < 1e-9
        % no uncertainty (point predictions instead of prediction set)
        c_y = center(Y);
        if max(c_y) <= c_y(ii) + 1e-9
            % class ii is predicted
            num = num +1;
            if any(ii == class_true,'all')
                numCorrect = numCorrect +1;
            end
        end

    else
        % check prediction set
        T_i = repmat((1:n_y)==ii, n_y, 1) - eye(n_y);
        problem.f = zeros(1,eta);
        problem.Aineq = -T_i*generators(Y);
        problem.bineq = T_i*center(Y);
        problem.lb = -1*ones(eta,1);
        problem.ub = ones(eta,1);
        problem.options = optionsOptim;
        [~,~,exitflag,~] = CORAlinprog(problem);
        if exitflag >0
            % class ii lies in the prediction set
            num = num +1;
            if any(ii == class_true,'all')
                % class ii is correct class
                numCorrect = numCorrect +1;
            end
        end
    end
end
end

function [num, numCorrect] = aux_computeNumClassesCP(classes_pred,class_true)
% compute the number of classes which are feasible for a given output set

num = length(classes_pred);
numCorrect = 0;
for class_pred = classes_pred
    if any(class_pred == class_true,'all')
        numCorrect = numCorrect +1;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
