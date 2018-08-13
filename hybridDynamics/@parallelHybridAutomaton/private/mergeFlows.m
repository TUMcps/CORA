function  res = mergeFlows(obj, flowList)
% mergeFlows - Merge the continious dynamics of several subcomponents to 
%              obtain the continous dynamic for the overall system
%
% Syntax:  
%    res = mergeFlows(obj, flowList)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    flowList - continous dynamics object for each subcomponent
%
% Outputs:
%    res - constructed continous dynamics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018 (NK, output instead of state for input binds)
% Last revision: ---

%------------- BEGIN CODE --------------

    numComps = length(flowList);

    % check user input
    isLinSys = false(numComps);
    for i = 1:numComps
        isLinSys(i) = isa(flowList{i},'linearSys');
    end

    if ~all(isLinSys)
        error("Only LINEAR parallel hybrid automata are currently supported.");
    end

    % Once multiple dynamics are supported, the following code should be
    % moved into its own method, perhaps as part of the 'linearSys' class.

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
       c = flow.c;

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

       % merge names with ampersands
       name = flow.name;
       if i == 1
           nameMerged = name;
       else
           nameMerged = [nameMerged,' & ',name];
       end
    end

    % construct resulting continious dynamics object
    res = linearSys(nameMerged,Amerged,Bmerged,cMerged);

end

%------------- END OF CODE --------------