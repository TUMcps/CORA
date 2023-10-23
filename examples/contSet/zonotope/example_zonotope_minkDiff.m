function completed = example_zonotope_minkDiff()
% example_zonotope_minkDiff - evaluates the performance of the Minkowski
% difference compared to polytopes
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
runs = 10; % number of repetitions for averaging

%% length of generators for the minuend to avoid too many empty Minkowski differences
orderFactor = subtrahendOrder/minuendOrder; % factor between 
lengthFactor = 10*orderFactor; %lengthFactor = averageMinuendGenLength/averageSubtrahendGenLength

%% initialize sets
Zm = cell(runs,1); % Minuend represented as zonotope
Zs = cell(runs,1); % Subtrahend represented as zonotope
Zres = cell(runs,4); % Minkowski difference represented as zonotope
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

% set considered types
%typeSet = {'under', 'over', 'approx', 'overCoarse', 'conZonotope', 'RaghuramanKoeln'};
%typeSet = {'under', 'over', 'approx', 'overCoarse', 'RaghuramanKoeln', 'RaghuramanKoeln_weighted', 'underIncremental'};
%typeSet = {'under', 'RaghuramanKoeln', 'RaghuramanKoeln_weighted', 'underIncremental_basic', 'underIncremental'};
%typeSet = {'under', 'over', 'overCoarse', 'dummy'};
typeSet = {'inner'};

% set number of considered types
types = length(typeSet);

%% compute Minkowski difference
for iType = 1:types

    % set type
    type = typeSet{iType};
    
    tStart = tic;
    %profile on
    for i = 1:length(Zm) 
        Zres{i,iType} = minkDiff(Zm{i},Zs{i},type);
    end
    %profile off
    %profile viewer
    t_zono = toc(tStart);
    t_average(iType) = t_zono/runs;
end
t_average


%% analyze results
emptySet = zeros(runs,length(Zres(1,:)));
order= zeros(runs,length(Zres(1,:)));
reducedOrder = zeros(runs,length(Zres(1,:)));
% loop over all computation types
for iType = 1:types
    % loop over all zonotopes
    for i = 1:length(Zres(:,iType)) 
        % is zonotope empty?
        if representsa(Zres{i,iType},'emptySet',1e-12)
            emptySet(i,iType) = 1;
        else
            % remove zero generators
            Zred = compact(Zres{i,iType},'zeros');
            G = generators(Zred);
            order(i,iType) = length(G(1,:))/n;
            if order(i,iType) ~= minuendOrder
                reducedOrder(i,iType) = 1;
            end
        end
    end
end

emptySets = sum(emptySet,1);
averageOrder = sum(order,1)./(runs-emptySets);

%% compute Minkowski difference for polytopes
tStart = tic;
for i = 1:length(Pm)
    Pres{i} = minkDiff(Pm{i}, Ps{i});
end
t_poly = toc(tStart);
t_averagePoly = t_poly/runs

%% analyze result
emptySetPoly = zeros(1,runs);
reducedOrderPoly = zeros(1,runs);
%for i = 1:length(Zres) 
for i = 1:0 
    % is polytope empty?
    if representsa(Pres{i},'emptySet')
        emptySetPoly(i) = 1;
    else
        %check size of new constraints
        Pred = removeRedundancies(Pres{i});
        P_h = halfspace(Pred);
        K = get(P_h,'K');
        constraintsNew = length(K);
        
        %check size of old constraints
        P_h = halfspace(Pm{i});
        K = get(P_h,'K');
        constraintsOld = length(K);
        
        %check size of constraints
        if constraintsNew < constraintsOld
            reducedOrderPoly(i) = 1;
        end
    end
end


% show number of empty sets
emptySets
if emptySet(1) ~= emptySetPoly
    disp('emptyCounter inconsistent')
end

% show reduced order
averageOrder
if any(reducedOrder(:,1) > reducedOrderPoly') % reduced number of halspace doe not necessarily reduce number of generators
    disp('reduced order inconsistent')
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
% init index of non-empty results
ind = 0; 
% loop over all zonotopes
for i = 1:length(Zres(:,iType)) 
    % polytope empty?
    if ~representsa(Zres{i,1},'emptySet')
        % loop over all computation types
        for iType = 1:types
            % sample containment of zonotope
            sc_Z(iType) = volume(Zres{i,iType});
        end
        if sc_Z(1) > 0 % under-approximation is reference
            % increment counter of nenempty sets
            ind = ind + 1
            % relative volume
            vol_rel(ind,:) = sc_Z/sc_Z(1);
            % normalized relative volume
            vol_relNorm(ind,:) = sc_Z.^(1/n)/sc_Z(1).^(1/n);
        end
    end
end

% average volumes
vol_rel_avg = sum(vol_rel,1)/ind
vol_relNorm_avg = sum(vol_relNorm,1)/ind


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
tStart = tic;
for i = 1:length(Zres) 
    %compute Minkowski difference
    if representsa(Zres{i},'emptySet')
        Zres_add{i} = Zs{i};
    else
        Zres_add{i} = Zres{i} + Zs{i};
    end
end
t_zono_add = toc(tStart);
t_average_add = t_zono_add/runs

%% Minkowski addition for polytopes
tStart = tic;
for i = 1:length(Pres) 
    Pres_add{i} = Pres{i} + Ps{i};
end
t_poly_add = toc(tStart)
t_PolyAverage_add = t_poly_add/runs

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
