function example_nonlinHybridDiscrete_reach_ARCH18_pancreas()


    % Initialization ------------------------------------------------------

    % Initial set
    int = interval([0;72.429;141.149;162.449;229.824;3.19;5.49;100.249;100.249;120;0;50], ...
                   [0;72.431;141.151;162.451;268.128;3.21;5.51;100.251;100.251;160;0;90]);

    R0 = zonotope(int);

    options.x0=center(R0); % initial state for simulation
    options.R0=R0; % initial state for reachability analysis

    % additional options
    options.startLoc = 3; % initial location
    options.finalLoc = -inf; % no final location
    options.tStart=0; % start time
    options.tFinal=720; % final time
%     options.tFinal = 200;
    options.intermediateOrder = 10;
    options.originContained = 0;
    options.timeStepLoc = repmat({0.5},[30,1]);

    options.zonotopeOrder=10; % zonotope order
    options.polytopeOrder=3; % polytope order
    options.reductionTechnique = 'girard';
    options.isHyperplaneMap=0;
    options.guardIntersect = 'zonoGirard';     % use constrained zonotopes to calculate guard intersections
    options.intersectInvariant = 0; 
    options.enclosureEnables = [1,2]; % choose enclosure method(s)
    options.filterLength = [5,7];

    options.taylorTerms=20;             % number of taylor terms for reachable sets
    options.timeStepAll=5;
    options.advancedLinErrorComp = 0;
    options.tensorOrder = 2;
    options.errorOrder = 10;
    options.reductionInterval = inf;
    options.loc=3;
    options.reductionTechnique='girard';
    options.maxError = 150e5*ones(12,1);
    options.oldError = zeros(12,1);
    
    options.plotType = {'b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y','b','m','g','r','y'};

    % specify hybrid automata
    HA = pancreas2_hp(); % automatically converted from SpaceEx
    locList = get(HA,'location');

    % divide guards in time-guards and state-guards
    stateGuard = cell(length(locList),1);
    timeGuard = cell(length(locList),1);
    
    for i = 1:length(locList)
    
        trans = get(locList{i},'transition');

        counter = 1;
        
        for k = 1:length(trans)
            guard = get(trans{k},'guard');
            ind = find(abs(guard.c)>0);
            if ind == 10
               stateGuard{i}{counter}.guard = guard;
               stateGuard{i}{counter}.target = get(trans{k},'target');
               counter = counter + 1;
            elseif ind == 11
               timeGuard{i}.guard = guard;
               timeGuard{i}.target = get(trans{k},'target');
            end
        end
    end
    
    
    
    
    % Simulation ----------------------------------------------------------
    
    % draw random points from the initial set as starting points for the
    % simulation
    fractionVertices = 0.5;
    N = 100;
    
    points = cell(N,1);
    locIDpoint = cell(N,1);
    counter = 1;
    
    for i = 1:length(points)
        if counter < fractionVertices * N
            points{i} = randPointExtreme(R0);
        else
            points{i} = randPoint(R0);
        end
        locIDpoint{i} = options.startLoc;
        counter = counter + 1;
    end
    
    % initialization
    t = options.tStart : options.timeStepAll : options.tFinal;
    dataSim = cell(length(t),1);
    
    % loop over all time periods
    for i = 1:length(t)   
        
        dataSim{i}.locs = locIDpoint;
        dataSim{i}.traj = cell(length(locIDpoint),1);
        
        points_ = cell(length(points),1);
        locIDpoint_ = cell(length(points),1);
        
        % loop over all points
        for j = 1:length(points)
            
            % set options
            sys = get(locList{locIDpoint{j}},'contDynamics');
            options.R0 = zonotope(points{j});
            options.tFinal = options.timeStepAll;
            options.tStart = 0;
            options.timeStep = options.timeStepLoc{locIDpoint{j}};
            options.U = zonotope(0);
            options.uTrans = 0;
            options.u = 0;
            
            % reachability analysis
            [~,simT,simX] = simulate(sys,options,options.tStart,options.tFinal,points{j},options);
            P = simX(end,:)';
            
            % calculate guard intersection (state variables)
            guard = stateGuard{locIDpoint{j}};
            setFinished = 0;
            
            for k = 1:length(guard)
                
               g = guard{k}.guard;
                   
               % point left current invariant
               if g.c' * P < g.d
                   
                   points_{j} = P;
                   locIDpoint_{j} = guard{k}.target;
                   setFinished = 1;
                   break;
                   
               end
            end
            
            % if the set stayed in the invariant
            if ~setFinished
               points_{j} = P;
               locIDpoint_{j} = locIDpoint{j};
            end
      
            % store data
            dataSim{i}.traj{j}.t = simT;
            dataSim{i}.traj{j}.x = simX;
        end
        
        % intersection with time guards
        for j = 1:length(points_)
            
            guard = timeGuard{locIDpoint_{j}};
            
            if i < length(t) && ~isempty(guard) && t(i+1) == -guard.guard.d
                locIDpoint_{j} = guard.target;
            end     
        end
        
        % update values
        points = points_;
        locIDpoint = locIDpoint_;
    end
    
    
    
    
    
    
    % Reachability Analysis -----------------------------------------------

    sets{1} = R0;
    locIDs{1} = options.startLoc;

    data = cell(length(t),1);
    
    % loop over all time periods
    for i = 1:length(t)
        
        sets_ = cell(length(sets),1);
        locIDs_ = cell(length(sets),1);
        
        data{i}.locs = locIDs;
        data{i}.sets = cell(length(sets),1);
        
        counter = 1;
        
        % loop over all parallel sets
        for j = 1:length(sets)
            
            % set options
            sys = get(locList{locIDs{j}},'contDynamics');
            options.R0 = sets{j};
            options.tFinal = options.timeStepAll;
            options.tStart = 0;
            options.timeStep = options.timeStepLoc{locIDs{j}};
            options.U = zonotope(0);
            options.uTrans = 0;
            
            % reachability analysis
            [Rcont,Rtp] = reach(sys,options); 
            R = Rtp{end}{1}.set;
            
            % calculate guard intersection (state variables)
            guard = stateGuard{locIDs{j}};
            setFinished = 0;
            
            for k = 1:length(guard)
                
               g = guard{k}.guard;
               
               % set intersects two invariants -> split
               if isIntersectingApprox(g,R)
                   
                   temp = intervalOverApprox(R);
                   Z1 = temp;
                   Z2 = temp;
                   
                   if g.c(10) > 0
                      Z2(10) = interval(g.d,supremum(Z2(10)));
                      Z1(10) = interval(infimum(Z1(10)),g.d);
                   else
                      Z2(10) = interval(infimum(Z2(10)),-g.d);
                      Z1(10) = interval(-g.d,supremum(Z1(10)));
                   end
                   
                   if isa(R,'zonotopeBundle')
                       
                       % set outside the invariant
                       Z = R.Z;
                       Z{end+1} = zonotope(Z1);
                       sets_{counter} = zonotopeBundle(Z);
                       locIDs_{counter} = guard{k}.target;
                       counter = counter + 1;
                       
                       % set inside the current invariant
                       Z = R.Z;
                       Z{end+1} = zonotope(Z2);
                       sets_{counter} = zonotopeBundle(Z);
                       locIDs_{counter} = locIDs{j};
                       counter = counter + 1;
                       
                   else
                       
                       % set outside the invariant
                       sets_{counter} = zonotopeBundle({R,zonotope(Z1)});
                       locIDs_{counter} = guard{k}.target;
                       counter = counter + 1;
                       
                       % set inside the current invariant
                       sets_{counter} = zonotopeBundle({R,zonotope(Z2)});
                       locIDs_{counter} = locIDs{j};
                       counter = counter + 1;
                       
                   end
                   
                   setFinished = 1;
                   
                   
               % set fully outside invariant    
               elseif g.c' * center(R) < g.d
                   
                   sets_{counter} = R;
                   locIDs_{counter} = guard{k}.target;
                   counter = counter + 1;
                   
                   setFinished = 1;
          
               end
            end
            
            % if the set stayed in the invariant
            if ~setFinished
               sets_{counter} = R;
               locIDs_{counter} = locIDs{j};
               counter = counter + 1;
            end
      
            % store data
            data{i}.sets{j} = Rcont;
            data{i}.loc{j} = locIDs{j};
        end
        
        
        % intersection with time guards
        for j = 1:length(sets_)
            
            guard = timeGuard{locIDs_{j}};
            
            if i < length(t) && ~isempty(guard) && t(i+1) == -guard.guard.d
                locIDs_{j} = guard.target;
            end     
        end
        
        % merge reachable sets
        temp = cell2mat(locIDs_);
        locAll = unique(temp);
        sets__ = cell(length(locAll),1);
        locIDs__ = cell(length(locAll),1);
        
        for j = 1:length(locAll)
            
            indices = find(temp == locAll(j));
            
            % more than one set in the same location -> merge
            if length(indices) > 1
                
                sets__{j} = mergeSets(sets_(indices),options);
                locIDs__{j} = locIDs_{indices(1)};
                
%                 inter = interval(sets_{indices(1)});
%                 sup = supremum(inter);
%                 infi = infimum(inter);
% 
%                 for k = 2:length(indices)
%                     inter = interval(sets_{indices(k)});
% 
%                     sup = max(supremum(inter),sup);
%                     infi = min(supremum(inter),infi);
%                 end
% 
%                 sets__{j} = zonotope(interval(infi,sup));
%                 locIDs__{j} = locIDs_{indices(1)};
                
            else
                sets__{j} = sets_{indices(1)};
                locIDs__{j} = locIDs_{indices(1)};
            end
        end
        
        % update values
        sets = sets__;
        locIDs = locIDs__;
        
        
        % Debug
        hold on
        
        for j = 1:length(data{i}.sets)
            tempSet = data{i}.sets{j};
            for k = 1:length(tempSet)
                temp = reduce(project(tempSet{k}{1},[11,10]),'girard',3);
                plot(temp,[1,2],options.plotType{data{i}.locs{j}});
            end
        end
        
        for j = 1:length(dataSim{i}.traj)
           plot(dataSim{i}.traj{j}.x(:,11), dataSim{i}.traj{j}.x(:,10),'k');
        end
        
    end


    % Visualization -------------------------------------------------------








end



% Auxiliary Functions -----------------------------------------------------

function R = mergeSets(Rall,options)

    % determine suitalbe othoganal basis
    D = determineBasis(Rall,options);
    
    Rtemp = cell(length(D),1);
    
    % loop over all different basis
    for i = 1:length(D)
       
        infi = inf * ones(size(D{i},1),1);
        sup = -inf * ones(size(D{i},1),1);

        % loop over all sets
        for j = 1:length(Rall)
               
            % transform the reachable set to the space defined by the
            % current basis and over-approximate the set with an interval
            temp = intervalOverApprox(D{i}'*Rall{j});
                   
            % update bounds
            if ~isempty(temp)
               infi = min(infi,infimum(temp));
               sup = max(sup,supremum(temp));
            end
        end

        % transform bounding box back to original space
        temp = interval(infi,sup);
        Rtemp{i} = D{i}*zonotope(temp);       
    end
    
    % construct final set
    if length(Rtemp) == 1
       R = Rtemp{1}; 
    else
       R = zonotopeBundle(Rtemp); 
    end
end


function D = determineBasis(Rall,options)

    D = cell(length(options.enclosureEnables),1);
    
    % loop over all enclosure methods
    for i = 1:length(D)
        
        switch options.enclosureEnables(i)
           
            % axis aligned bounding box
            case 1
                
                if isa(Rall{1},'zonotopeBundle')
                    D{i} = eye(size(Rall{1}.Z{1}.Z,1));
                else
                    D{i} = eye(size(Rall{1}.Z,1));
                end
                
                
                
            % Principal Component Analysis (generators)
            case 2               
                
                % concatenate all generators
                G = [];
                
                for j = 1:length(Rall)
                   if isa(Rall{j},'zonotopeBundle')
                      for k = 1:Rall{j}.parallelSets
                         G = [G, Rall{j}.Z{k}.Z(:,2:end)]; 
                      end
                   else
                       G = [G, Rall{j}.Z(:,2:end)];
                   end
                end                
                
                % calcualte an orthogonal transformation matrix using PCA
                [D{i},~,~] = svd(G); 

        end               
    end
end


function res = intervalOverApprox(R)

    if isa(R,'zonotopeBundle')
       
       temp = interval(R.Z{1});
       sup = supremum(temp);
       infi = infimum(temp);
        
       for i = 2:R.parallelSets
           temp = interval(R.Z{i});
           sup = min(sup,supremum(temp));
           infi = max(infi,infimum(temp));
       end
       
       for j = 1:length(sup)
          if abs(sup(j) - infi(j)) < 1e-12
             val = (sup(j) + infi(j)) / 2;
             sup(j) = val;
             infi(j) = val;
          end
       end
       
       res = interval(infi,sup);
    else
       res = interval(R); 
    end

end