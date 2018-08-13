function mergedBC = automatonProduct(BC1,BC2)
% Compute the automaton product of two Base Component Instances
% (The input BCs should be brought to the same variable space first)
% RETURNS: a single Base Component Instance

% Meta-information are combined to stay informative & readable
mergedBC.id = BC1.id + " X " + BC2.id;
mergedBC.name = BC1.name + " X " + BC2.name;

% Commence main task: building the cross-product of the state sets.
S1 = BC1.States;
numS1 = numel(S1);
S2 = BC2.States;
numS2 = numel(S2);

% Build cross product S1 X S2
S1XS2 = struct('id',cell(1,numS1*numS2));

for i = 1:numS1
    for j = 1:numS2
        % combined State (S1(i),S2(j)) is assigned to S1XS2((i-1)*|S2| + j)
        idx = (i-1)*numS2 + j;
        
        % Combine Meta-Information
        S1XS2(idx).id = S1(i).id + " X " + S2(j).id;
        S1XS2(idx).name = S1(i).name + " X " + S2(j).name;

        % Combine Flow, by concatenating the parsed flow equations.
        % (reminder: flow expressions are stored vertically)
        vn = [S1(i).Flow.varNames; S2(j).Flow.varNames];
        exprs = [S1(i).Flow.expressions; S2(j).Flow.expressions];
        % combine equation Text, to keep it sensible and readable
        text = S1(i).Flow.Text +newline+"&&"+newline+ S2(j).Flow.Text;
        
        S1XS2(idx).Flow = struct('varNames',vn,'expressions',exprs,'Text',text);
        
        
        % Combine Invariant, by building the cross product of its polytope.
        % This can be achieved, by concatenating the parsed expressions.
        % (reminder: polytope expressions are stored vertically)
        ineq = [S1(i).Invariant.inequalities; S2(j).Invariant.inequalities];
        eq = [S1(i).Invariant.equalities; S2(j).Invariant.equalities];
        % combine equation Text, to keep it sensible and readable
        text = S1(i).Invariant.Text +newline+"&&"+newline+ S2(j).Invariant.Text;
        
        S1XS2(idx).Invariant = struct('inequalities',ineq,'equalities',eq,'Text',text);
        
        % Combine Transition Sets of S1(i) & S2(j)
        T1 = S1(i).Trans;
        T2 = S2(j).Trans;
        % (update information that changes, then concatenate sets)
        
        % Foreach (S1_i -> S1_x) in T1, T1XT2 has ((S1_i,S2_j) -> (S1_x,S2_j))
        for t = 1:numel(T1)
            x = T1(t).destination;
            % compute index of new destination (S1_x,S2_j)
            T1(t).destination = (x-1)*numS2 + j;
            
            % Guard and Reset function can be left identical.
            % Polytopes naturally expand if new dimensions are added.
            % Reset for unmentioned variables is defined implicitly.
        end
        
        % analogous for (S2_j -> S2_x) T2
        for t = 1:numel(T2)
            x = T2(t).destination;
            % compute index of new destination (S1_i,S2_x)
            T2(t).destination = (i-1)*numS2 + x;
            
            % guard and reset function can be left identical,
            % since polytopes naturally expand if new dimensions are added
            % and reset functions are implicitly defined for unmentioned variables
        end
        
        % Concatenate tweaked Transition sets
        S1XS2(idx).Trans = [T1 T2];
        
    end
end

% add new state set to output component
mergedBC.States = S1XS2;

end