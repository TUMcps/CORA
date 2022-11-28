function sys = mergeFlows(pHA,flowList,loc)
% mergeFlows - merges the continuous dynamics of several subcomponents to 
%    obtain the continous dynamics for the overall system
%
% Syntax:  
%    sys = mergeFlows(pHA,flowList,loc)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    flowList - continous dynamics object for each subcomponent
%    loc - indices of the current location
%
% Outputs:
%    sys - constructed continous dynamics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schoepfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018 (NK, output instead of state for input binds)
% Last revision: ---

%------------- BEGIN CODE --------------

    % number of components
    numComps = length(flowList);

    % check user input
    isLinSys = false(numComps,1);
    isNonlinSys = false(numComps,1);
    for i=1:numComps
        isLinSys(i) = isa(flowList{i},'linearSys');
        isNonlinSys(i) = isa(flowList{i},'nonlinearSys');
    end

    % merge flows according to the dynamics
    if all(isLinSys)
        sys = mergeFlowsLinearSys(pHA,flowList);
    elseif all(isNonlinSys)
        sys = mergeFlowsNonlinearSys(pHA,flowList,loc);
    else
        throw(CORAerror('CORA:specialError',...
            ['Only "linearSys" and "nonlinearSys" objects are currently ', ...
            'supported for parallel hybrid automata are currently supported!']));
    end 
end


% Auxiliary Functions -----------------------------------------------------

function res = mergeFlowsLinearSys(pHA,flowList)

    numComps = length(flowList);

    % allocate merged dynamics
    Amerged = zeros(pHA.numStates,pHA.numStates);
    Bmerged = zeros(pHA.numStates,pHA.numInputs);
    cMerged = zeros(pHA.numStates,1);
    name = cell(numComps,1);

    % loop over all subcomponents
    for i = 1:numComps

       % get object properties
       flow = flowList{i};
       stateBinds = pHA.bindsStates{i};
       inputBinds = pHA.bindsInputs{i};

       name{i,1} = flow.name;
       A = flow.A;
       B = flow.B;
       if isempty(flow.c)
           c = zeros(size(A,1),1);
       else
           c = flow.c;
       end

       % constant input vector c
       cMerged(stateBinds) = cMerged(stateBinds) + c;

       % system matrix A
       Amerged(stateBinds,stateBinds) = A;

       % input matrix B
       for j = 1:size(inputBinds,1)

          if inputBinds(j,1) == 0         % global input
              Bmerged(stateBinds,inputBinds(j,2)) = ...
                  Bmerged(stateBinds,inputBinds(j,2)) + B(:,j);

          else                            % input = output of other component

              % equation y = C*x + D*u + k
              tempFlow = flowList{inputBinds(j,1)};
              tempStateBinds = pHA.bindsStates{inputBinds(j,1)};
              tempInputBinds = pHA.bindsInputs{inputBinds(j,1)};

              % if a system has no matrices C/D/k, we assume the states are
              % given as outputs for parallelization (just as for nonlinear
              % systems)
              if isscalar(tempFlow.C) && tempFlow.C == 1 ...
                      && isempty(tempFlow.D) && isempty(tempFlow.k)
                  C = eye(tempFlow.dim);
              else 
                  C = tempFlow.C;
              end
              
              
              % part with matrix C
              Amerged(stateBinds,tempStateBinds) = ...
                  Amerged(stateBinds,tempStateBinds) + B(:,j)*C(inputBinds(j,2),:);
              
              % part with offset vector k
              if ~isempty(tempFlow.k)
                  cMerged(stateBinds) = cMerged(stateBinds) + B(:,j)*tempFlow.k(inputBinds(j,2));
              end
              
              % part with throughput matrix D
              if ~isempty(tempFlow.D)
                  D = tempFlow.D;
                  d = D(inputBinds(j,2),:);
                  ind1 = find(tempInputBinds(:,1) == 0);
                  ind2 = setdiff(1:size(tempInputBinds,1),ind1);

                  % check if the D matrix is valid
                  if any(d(ind2))
                      throw(CORAerror('CORA:specialError',...
                          ['It is not allowed for the throughput matrix D '...
                             'to point to inputs that are defined by the '...
                             'output of other subsystems, since this would '...
                             'otherwise lead to infinite loops!']));
                  end
                  
                  % construct the merged B matrix from the throughput
                  Bmerged(stateBinds,tempInputBinds(ind1,2)) = ...
                        Bmerged(stateBinds,tempInputBinds(ind1,2)) + ...
                        B(:,j)*d(ind1);
              end
          end
       end
    end

    % allocate merged name (remove all default names '')
    namemerged = strjoin(name(cellfun(@(x)~isempty(x),name,'UniformOutput',true)),' x ');

    % construct resulting continuous dynamics object
    res = linearSys(namemerged,Amerged,Bmerged,cMerged);
end

function res = mergeFlowsNonlinearSys(pHA,flowList,loc)

    numComps = length(flowList);

    % construct symbolic state vector and input vector
    x = sym('x',[pHA.numStates,1]);
    u = sym('u',[pHA.numInputs,1]);
    
    % initialize dynamic function
    f = sym(zeros(pHA.numStates,1));
    
    % loop over all subcomponents
    for i = 1:numComps

       % get object properties
       flow = flowList{i};
       stateBinds = pHA.bindsStates{i};
       inputBinds = pHA.bindsInputs{i};

       % construct input vector for this subcomponent
       u_ = sym(zeros(size(inputBinds,1),1));
       
       for j = 1:size(inputBinds,1)
          
           if inputBinds(j,1) == 0    % global input
               
               u_(j) = u(inputBinds(j,2));

           else                       % input = output of other component
           
               tempStateBinds = pHA.bindsStates{inputBinds(j,1)};
               u_(j) = x(tempStateBinds(inputBinds(j,2)));
               
           end
       end
       
       % construct flow function for this subcomponent
       f(stateBinds,:) = flow.mFile(x(stateBinds),u_);
    end

    % construct resulting continuous dynamics object
    name = ['location',strrep(strcat(num2str(loc')),' ','')];
    funHan = matlabFunction(f,'Vars',{x,u});
    
    res = nonlinearSys(name,funHan,pHA.numStates,pHA.numInputs);
end

%------------- END OF CODE --------------