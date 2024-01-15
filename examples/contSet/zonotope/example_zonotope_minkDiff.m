function completed = example_zonotope_minkDiff()
% example_zonotope_minkDiff - evaluates the performance of the Minkowski
%    difference compared to polytopes
%
% Syntax:
%    completed = example_zonotope_minkDiff()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-August-2015 
% Last update:   30-June-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Set parameters to randomly generate zonotopes
n = 6; % dimension
minuendOrder = 2; % order of minuend
subtrahendOrder = 2; % order of subtrahend
runs = 2; % number of repetitions for averaging

%% length of generators for the minuend to avoid too many empty Minkowski differences
orderFactor = subtrahendOrder/minuendOrder; % factor between 
lengthFactor = 10*orderFactor; %lengthFactor = averageMinuendGenLength/averageSubtrahendGenLength

%% set considered types
%typeSet = {'inner', 'outer', 'approx', 'overCoarse', 'conZonotope', 'RaghuramanKoeln'};
%typeSet = {'inner', 'outer', 'approx', 'overCoarse', 'RaghuramanKoeln', 'RaghuramanKoeln_weighted', 'underIncremental'};
%typeSet = {'inner', 'RaghuramanKoeln', 'RaghuramanKoeln_weighted', 'underIncremental_basic', 'underIncremental'};
%typeSet = {'inner', 'outer', 'overCoarse', 'dummy'};
typeSet = {'inner'};

% set number of considered types
numTypes = length(typeSet);

%% initialize sets
Zm = cell(runs,1); % Minuend represented as zonotope
Zs = cell(runs,1); % Subtrahend represented as zonotope
Zres = cell(runs,numTypes); % Minkowski difference represented as zonotope
Pm = cell(runs,1); % Minuend represented as polytope
Ps = cell(runs,1); % Subtrahend represented as polytope
Pres = cell(runs,1); % Minkowski difference represented as polytope

%% generate zonotopes and polytopes
for i = 1:runs
    %generate random minuend
    Z_tmp = zonotope.generateRandom('Dimension',n,'NrGenerators',minuendOrder*n);
    Zm{i} = enlarge(Z_tmp,lengthFactor);
    %generate random subtrahend
    Zs{i} = zonotope.generateRandom('Dimension',n,'NrGenerators',subtrahendOrder*n);
    
    %convert zonotopes to polytopes
    Pm{i} = polytope(Zm{i});
    Ps{i} = polytope(Zs{i});
end

%% compute Minkowski difference for zonotopes
t_average = zeros(numTypes,1);
for iType = 1:numTypes

    % select type
    type = typeSet{iType};
    
    tStart = tic;
    %profile on
    for i = 1:runs
        Zres{i,iType} = minkDiff(Zm{i},Zs{i},type);
    end
    %profile off
    %profile viewer
    t_zono = toc(tStart);
    t_average(iType) = t_zono/runs;
    disp("Average time of Mink.diff. for zonotopes (type = " + type + ...
        "): " + t_average(iType));
end

%% analyze results
isEmptySet = false(runs,numTypes);
order = zeros(runs,numTypes);
isReducedOrder = false(runs,numTypes);
% loop over all computation types
for iType = 1:numTypes
    % loop over all zonotopes
    for i = 1:runs
        % is zonotope empty?
        if representsa_(Zres{i,iType},'emptySet',1e-12)
            isEmptySet(i,iType) = true;
        else
            % remove zero generators
            Zred = compact(Zres{i,iType},'zeros');
            order(i,iType) = size(generators(Zred),2)/n;
            isReducedOrder(i,iType) = order(i,iType) ~= minuendOrder;
        end
    end
end
numEmptySets = sum(isEmptySet,1);
averageOrder = sum(order,1)./(runs-numEmptySets);

%% compute Minkowski difference for polytopes
tStart = tic;
for i = 1:length(Pm)
    disp(i)
    Pres{i} = minkDiff(Pm{i}, Ps{i});
end
t_poly = toc(tStart);
t_averagePoly = t_poly/runs;
disp("Average time of Mink.diff. for polytopes: " + t_averagePoly);

%% analyze result
isEmptySetPoly = false(1,runs);
reducedOrderPoly = false(1,runs);
for i = 1:runs
    % is polytope empty?
    if representsa_(Pres{i},'emptySet',1e-12)
        isEmptySetPoly(i) = true;
    else
        %check size of new constraints
        Pred = compact(Pres{i});
        constraintsNew = length(Pred.b) + length(Pred.be);
        
        %check size of old constraints
        constraintsOld = length(Pm{i}.b) + length(Pm{i}.be);
        
        %check size of constraints
        reducedOrderPoly(i) = constraintsNew < constraintsOld;
    end
end


% show number of empty sets
disp("Number of empty results: " + numEmptySets);
if isEmptySet(1) ~= isEmptySetPoly
    disp('emptyCounter inconsistent');
end

% show reduced order (reduced number of halfspaces does not necessarily
% reduce number of generators)
disp("Average order of results: " + averageOrder);
if any(isReducedOrder(:,1) > reducedOrderPoly')
    disp('reduced order inconsistent');
end


% %% compute relative volumes (exact -- only for small dimensions)
% % init index of non-empty results
% ind = 0; 
% % loop over all zonotopes
% for i = 1:length(Zres(:,iType)) 
%     % polytope empty?
%     if ~representsa(Zres{i,1},'emptySet')
%         % increment counter of nenempty sets
%         ind = ind + 1;
%         % compute volume of polytope
%         vol_P = volume(Pres{i});
%         % loop over all computation types
%         for iType = 1:3
%             % volume of zonotope
%             try
%                 vol_Z(iType) = volume(Zres{i,iType});
%             catch
%                 disp('volume could not be computed')
%             end
%         end
%         % relative volume
%         vol_rel(ind,:) = abs(vol_Z - vol_P)/vol_P;
%         % normalized relative volume
%         vol_relNorm(ind,:) = abs(vol_Z.^(1/dim) - vol_P.^(1/dim))/vol_P.^(1/dim);
%     end
% end
% 
% % average volumes
% vol_rel_avg = sum(vol_rel,1)/ind
% vol_relNorm_avg = sum(vol_relNorm,1)/ind

% %% compute relative volumes (interval approximation)
% % init index of non-empty results
% ind = 0; 
% % loop over all zonotopes
% for i = 1:length(Zres(:,iType)) 
%     % polytope empty?
%     if ~representsa(Zres{i,1},'emptySet')
%         % increment counter of nenempty sets
%         ind = ind + 1;
%         % compute volume of polytope
%         vol_P = volume(interval(Pres{i}));
%         % loop over all computation types
%         for iType = 1:3
%             % volume of zonotope
%             vol_Z(iType) = volume(interval(Zres{i,iType}));
%         end
%         % relative volume
%         vol_rel(ind,:) = abs(vol_Z - vol_P)/vol_P;
%         % normalized relative volume
%         vol_relNorm(ind,:) = abs(vol_Z.^(1/dim) - vol_P.^(1/dim))/vol_P.^(1/dim);
%     end
% end
% 
% % average volumes
% vol_rel_avg = sum(vol_rel,1)/ind
% vol_relNorm_avg = sum(vol_relNorm,1)/ind


% %% compute relative volumes (Monte Carlo approximation)
% % init index of non-empty results
% ind = 0; 
% % nr of Samples
% N = 1e3;
% % loop over all zonotopes
% for i = 1:length(Zres(:,iType)) 
%     % polytope empty?
%     if ~representsa(Zres{i,1},'emptySet')
%         % sample points from over-approximation
%         S = randPoint(polytope(Zres{i,4}),N);
%         % compute volume of polytope
%         sc_P = aux_sampleContainment(Pres{i},S);
%         % loop over all computation types
%         for iType = 1:types
%             % sample containment of zonotope
%             sc_Z(iType) = aux_sampleContainment(Zres{i,iType},S);
%         end
%         if sc_P > 0
%             % increment counter of nenempty sets
%             ind = ind + 1
%             % relative volume
%             vol_rel(ind,:) = sc_Z/sc_P;
%             % normalized relative volume
%             vol_relNorm(ind,:) = sc_Z.^(1/dim)/sc_P.^(1/dim);
%         end
%     end
% end

%% compute relative volumes (relative to under-approximation)

sc_Z = zeros(runs,numTypes);
vol_rel = zeros(runs,numTypes);
vol_relNorm = zeros(runs,numTypes);
% loop over all zonotopes
for i = 1:runs
    % loop over all computation types
    for iType = 1:numTypes
        % compute volume of zonotope
        sc_Z(i,iType) = volume(Zres{i,iType});

    end

    if sc_Z(i,1) > 0 % under-approximation is reference
        % relative volume
        vol_rel(i,:) = sc_Z(i,:)/sc_Z(i,1);
        % normalized relative volume
        vol_relNorm(i,:) = sc_Z(i,:).^(1/n)/sc_Z(1,iType).^(1/n);
    end

    % compute and display average volumes
    vol_rel_avg = mean(vol_rel(i,:));
    vol_relNorm_avg = mean(vol_relNorm(i,:));
    disp("Run " + i);
    disp("..average relative volume: " + vol_rel_avg);
    disp("..average relative volume (norm): " + vol_relNorm_avg);
end


% %visual inspection
% % if ~isempty(indices)
% %     for i = indices
% for i=1:100
%         Zres{i}
%         Pres{i}
%         if ~representsa(Zres{i},'emptySet')
%             figure
%             plot(Zres{i},[1 2],'r');
%             plot(Pres{i},[1 2],'b'); %red should not be seen
%             plot(Zm{i},[1 2],'g');
%             plot(Zres{i} + Zs{i},[1 2],'r');
%             plot(Pres{i} + Ps{i},[1 2],'b'); %red and maybe green should not be seen
%         end
%     end
% % end


%% Minkowski addition
Zres_add = cell(runs,1);
tStart = tic;
for i = 1:runs
    %compute Minkowski addition
    Zres_add{i} = Zres{i} + Zs{i};
end
t_zono_add = toc(tStart);
t_zonoAverage_add = t_zono_add/runs;
disp("Average time of Mink.add. for zonotopes: " + t_zonoAverage_add);

%% Minkowski addition for polytopes
Pres_add = cell(runs,1);
tStart = tic;
for i = 1:runs
    disp(i)
    Pres_add{i} = Pres{i} + Ps{i};
end
t_poly_add = toc(tStart);
t_polyAverage_add = t_poly_add/runs;
disp("Average time of Mink.add. for polytopes: " + t_polyAverage_add);


completed = true;

end


% Auxiliary functions -----------------------------------------------------

% computes the fraction of samples inside a set
function sc = aux_sampleContainment(set,S)
    % init res
    res = zeros(length(S(1,:)),1);
    % loop over all points
    for i = 1:length(S(1,:))
        res(i) = in(set,S(:,i));
    end
    sc = sum(res)/length(res);
end


% ------------------------------ END OF CODE ------------------------------
