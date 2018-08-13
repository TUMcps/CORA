function [obj] = reach(obj,options)
% reach - computes the reachable set for a parallel hybrid automaton
%
% Syntax:  
%    [obj] = reach(obj,options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    options - options for reachability analysis
%
% Outputs:
%    obj - parallel hybrid automaton object with stored reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      04-July-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % load data from options
    tFinal = options.tFinal;                    % final time
    finalLoc = cell2mat(options.finalLoc);      % final location

    % initialize lists for open sets
    listTstart{1} = options.tStart;    % intermediate time at transitions
    listLoc{1} = options.startLoc;     % actual location
    listBlockedLoc{1} = repmat({0},[1,length(options.startLoc)]);
    listR0{1} = options.R0;

    % determine which inputs are specified locally and which globally
    inputMap = parseInputMap(obj,options);


    % loop over the different locations 
    while ~isempty(listLoc) && listTstart{1} < tFinal && ...
          ~all(finalLoc == cell2mat(listLoc{1}))
        
        % select first entry from list
        loc = listLoc{1};
        tStart = listTstart{1};
        blockedLoc = listBlockedLoc{1};
        R0 = listR0{1};
        
        listLoc = listLoc(2:end);
        listTstart = listTstart(2:end);
        listBlockedLoc = listBlockedLoc(2:end);
        listR0 = listR0(2:end);

        % construct new location with local Automaton Product
        currLocObj = locationProduct(obj,loc);

        % potentially convert the invariant set and the guard set to a
        % different set representation
        currLocObj = convGuard(currLocObj,options);

        % construct the input set for the current location
        options.U = mergeInputSet(obj,inputMap,loc,options);
        options.u = mergeInputVector(obj,inputMap,loc,options.uGlob,options.uCompLoc);
        options.uTrans = mergeInputVector(obj,inputMap,loc,options.uGlobTrans,options.uCompLocTrans);
        
        %obtain factors for initial state and input solution
        r = options.timeStep;
        
        for i=1:(options.taylorTerms+1)
            options.factor(i)= r^(i)/factorial(i);    
        end
        
        % reachability analysis for this location
        [TP,R,nextLoc,blockedLoc,Rjump,Rcont] = reach(currLocObj,tStart,R0,blockedLoc,options);
        
        
        % add resuling new branches to list
        listLoc = [listLoc,transpose(nextLoc)];
        listBlockedLoc = [listBlockedLoc, transpose(blockedLoc)];
        listTstart = [listTstart,num2cell(TP.tMin')];
        listR0 = [listR0,Rjump];
        
        % store the results
        if isempty(obj.result) || ~isfield(obj.result,'reachSet')
            obj.result.reachSet{1}.R = Rcont;
            obj.result.reachSet{1}.location = loc;
            obj.result.reachSet{1}.time = TP;
            obj.result.reachSet{1}.Rint = R;         
        else
            obj.result.reachSet{end+1}.R = Rcont;
            obj.result.reachSet{end}.location = loc;
            obj.result.reachSet{end}.time = TP;
            obj.result.reachSet{end}.Rint = R;
        end     
    end
    
    
    
    
% Auxiliary Functions -----------------------------------------------------

function U = mergeInputSet(obj,inputMap,loc,options)

    % initialize input zonotope center and generators
    Z = zeros(obj.numInputs,1);
    
    % loop over all used components
    comp = unique(inputMap(:,1));
    
    for i = 1:length(comp)
        
       % find indizes of component in global input set
       ind = find(inputMap(:,1) == comp(i));
       
       % get input set for this component
       if comp(i) ~= 0
            Utemp = options.UCompLoc{comp(i)}{loc{comp(i)}};
       else
            Utemp = options.UGlob;
       end
       
       % add generators to the global input set
       Z(ind,1) = Utemp.Z(inputMap(ind,2),1);
       
       if size(Utemp.Z,2) > 1
           G = zeros(obj.numInputs,size(Utemp.Z,2)-1);
           G(ind,:) = Utemp.Z(inputMap(ind,2),2:end);

           Z = [Z,G];
       end
    end
    
    % construct final zonotope object
    U = zonotope(Z);

%------------- END OF CODE --------------