function  res = mergeFlows(obj, flowList, loc)
% mergeFlows - Merge the continious dynamics of several subcomponents to 
%              obtain the continous dynamic for the overall system
%
% Syntax:  
%    res = mergeFlows(obj, flowList)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    flowList - continous dynamics object for each subcomponent
%    loc - indices of the current location
%
% Outputs:
%    res - constructed continous dynamics object
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

    numComps = length(flowList);

    % check user input
    isLinSys = false(numComps,1);
    isNonlinSys = false(numComps,1);
    for i = 1:numComps
        isLinSys(i) = isa(flowList{i},'linearSys');
        isNonlinSys(i) = isa(flowList{i},'nonlinearSys');
    end

    % merge flows according to the dynamics
    if all(isLinSys)
        res = mergeFlowsLinearSys(obj, flowList);
    elseif all(isNonlinSys)
        res = mergeFlowsNonlinearSys(obj, flowList, loc);
    else
        error(['Only "linearSys" and "nonlinearSys" objects are currently ', ...
       'supported for parallel hybrid automata are currently supported!']);
    end 
end


% Auxiliary Functions -----------------------------------------------------

function res = mergeFlowsLinearSys(obj, flowList)

    numComps = length(flowList);

    % allocate merged dynamics
    Amerged = zeros(obj.numStates,obj.numStates);
    Bmerged = zeros(obj.numStates,obj.numInputs);
    cMerged = zeros(obj.numStates,1);

    % loop over all subcomponents
    for i = 1:numComps

       % get object properties
       flow = flowList{i};
       stateBinds = obj.bindsStates{i};
       inputBinds = obj.bindsInputs{i};

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
              tempStateBinds = obj.bindsStates{inputBinds(j,1)};
              tempInputBinds = obj.bindsInputs{inputBinds(j,1)};
              C = tempFlow.C;

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
                      error(['It is not allowed for the throughput matrix D '...
                             'to point to inputs that are defined by the '...
                             'output of other subsystems, since it would '...
                             'otherwise be able to construct infinite loops!']);
                  end
                  
                  % construct the merged B matrix from the throughput
                  Bmerged(stateBinds,tempInputBinds(ind1,2)) = ...
                        Bmerged(stateBinds,tempInputBinds(ind1,2)) + ...
                        B(:,j)*d(ind1);
              end
          end
       end
    end

    % construct resulting continious dynamics object
    res = linearSys(Amerged,Bmerged,cMerged);
end

function res = mergeFlowsNonlinearSys(obj, flowList, loc)

    numComps = length(flowList);

    % construct symbolic state vector and input vector
    x = sym('x',[obj.numStates,1]);
    u = sym('u',[obj.numInputs,1]);
    
    % initialize dynamic function
    f = sym(zeros(obj.numStates,1));
    
    % loop over all subcomponents
    for i = 1:numComps

       % get object properties
       flow = flowList{i};
       stateBinds = obj.bindsStates{i};
       inputBinds = obj.bindsInputs{i};

       % construct input vector for this subcomponent
       u_ = sym(zeros(size(inputBinds,1),1));
       
       for j = 1:size(inputBinds,1)
          
           if inputBinds(j,1) == 0    % global input
               
               u_(j) = u(inputBinds(j,2));

           else                       % input = output of other component
           
               tempStateBinds = obj.bindsStates{inputBinds(j,1)};
               u_(j) = x(tempStateBinds(inputBinds(j,2)));
               
           end
       end
       
       % construct flow function for this subcomponent
       f(stateBinds,:) = flow.mFile(x(stateBinds),u_);
    end

    % construct resulting continuous dynamics object
    name = ['location',strrep(strcat(num2str(loc')),' ','')];
    funHan = matlabFunction(f,'Vars',{x,u});
    
    res = nonlinearSys(name,funHan,obj.numStates,obj.numInputs);
end

%------------- END OF CODE --------------