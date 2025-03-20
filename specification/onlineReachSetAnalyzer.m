classdef onlineReachSetAnalyzer < matlab.mixin.Copyable
% onlineReachSetAnalyzer - data structure for incrementally verifying a given STL formula
%
% Syntax:
%    analyzer = onlineReachSetAnalyzer(phi,propagationFrequency)
%
% Inputs:
%    phi - STL formula to verify
%    propagationFrequency - number of observations to accumulate before propagating
%
% Outputs:
%    analyzer - generated onlineReachSetAnalyzer object
%
% Example:
%    x = stl('x',2);
%    eq = until(x(1) < 6,x(2) > 10,stlInterval(0,1,true,false));
%    analyzer = onlineReachSetAnalyzer(eq,20)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: onlineAnalyzerNode

% Authors:       Florian Lercher
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    propagationFrequency (1,1) {mustBeInteger,mustBeNonnegative} = 1;
    nodes (1,:) onlineAnalyzerNode = onlineAnalyzerNode.empty;
    leaves (1,:) {mustBeInteger,mustBeNonnegative} = [];
    numObservations (1,1) {mustBeInteger,mustBeNonnegative} = 0;

    propagated (1,:) logical = logical.empty;
end

methods
    % constructor
    function obj = onlineReachSetAnalyzer(phi,propagationFrequency,dropPropagated)
        if nargin < 3
            % control whether to drop the propagated part of signals
            % this should only be turned off for debugging and visualization
            dropPropagated = true;
        end
        if nargin < 2
            propagationFrequency = 1;
        end
        cache = containers.Map;
        phi = desugar(phi);
        parse(phi,0);
        obj.propagated = false(1,length(obj.nodes));
        obj.leaves = setdiff(1:length(obj.nodes),vertcat(obj.nodes.parents));
        obj.propagationFrequency = propagationFrequency;

        % helper function to parse the STL formula into a syntax DAG
        function parse(phi,parent)
            % check cache
            k = formattedDisplayText(phi);
            if cache.isKey(k)
                obj.nodes(cache(k)).addParent(parent);
                obj.nodes(parent).addChild(cache(k));
                return
            end

            % no node for phi --> create one and store in cache
            myPos = length(obj.nodes) + 1;
            % need this check for the first iteration where we pass in 0 as parent
            if parent ~= 0
                obj.nodes(parent).addChild(myPos);
            end
            cache(k) = myPos;
            switch phi.type
                % atom
                case {'true','false'}
                    obj.nodes(myPos) = onlineAnalyzerNode(phi.type,myPos,parent,dropPropagated);
                case 'variable'
                    if isa(phi.lhs,'atomicProposition')
                        obj.nodes(myPos) = onlineAnalyzerNode('AP',myPos,parent,dropPropagated,phi.lhs);
                    else
                        throw(CORAerror('CORA:wrongInputInConstructor','STL type variable, but lhs is not an AP'));
                    end    
                    % comparator
                case {'<','<=','>','>='}
                    ap = atomicProposition(convert2set(phi));
                    obj.nodes(myPos) = onlineAnalyzerNode('AP',myPos,parent,dropPropagated,ap);
                    % boolean
                case '~'
                    obj.nodes(myPos) = onlineAnalyzerNode(phi.type,myPos,parent,dropPropagated);
                    parse(phi.lhs,myPos);
                case {'&','|'}
                    obj.nodes(myPos) = onlineAnalyzerNode(phi.type,myPos,parent,dropPropagated);
                    parse(phi.lhs,myPos);
                    parse(phi.rhs,myPos);
                    % temporal
                case 'finally'
                    obj.nodes(myPos) = onlineAnalyzerNode(phi.type,myPos,parent,dropPropagated,phi.interval);
                    parse(phi.lhs,myPos);
                case 'until'
                    obj.nodes(myPos) = onlineAnalyzerNode(phi.type,myPos,parent,dropPropagated,phi.interval);
                    parse(phi.lhs,myPos);
                    parse(phi.rhs,myPos);
                otherwise
                    throw(CORAerror('CORA:wrongInputInConstructor',strcat('unknown STL type: ',phi.type)));
            end
        end
    end

    % observe a set in a given interval and location
    % makeValid is an interval indicating for which time the observation are final
    function observeSet(obj,set,interval,loc,makeValid)
        % if the interval is empty, we cannot do anything with this observation
        if isemptyobject(interval)
            return;
        end

        % convert set into zonotope if it is a zonotope bundle
        if isa(set,'zonoBundle')
            set = zonotope(set);
        end

        if nargin < 5
            tentative = false;
            makeValid = stlInterval();
        else
            tentative = true;
        end
        
        % add the new observation to the leaves
        for l = obj.leaves
            obj.nodes(l).observeSet(set,interval,loc,tentative);
        end

        % trigger the propagation if necessary
        obj.numObservations = obj.numObservations + 1;
        if mod(obj.numObservations,obj.propagationFrequency) == 0
            obj.propagate(makeValid);
        end
    end

    % trigger propagation regardless of the propagation frequency
    function forcePropagation(obj)
        obj.propagate(stlInterval(0,inf));
    end

    % get the current verdict on satisfaction of the top-level formula
    function verdict = getVerdict(obj)
        verdict = obj.nodes(1).signal.at(0);
    end

    % plot the graph structure
    function han = plot(obj)
        if length(obj.nodes) == 1
            % plot graph
            han = plot(graph(1,1,'omitselfloops'),'NodeLabel',{formattedDisplayText(obj.nodes{1})});
        else
            % convert to Matlab directed graph
            edgeEnd = vertcat(obj.nodes.children)';
            numChildren = cellfun(@length,{obj.nodes.children});
            edgeStart = repelem(1:length(obj.nodes),numChildren);
            G = digraph(edgeStart,edgeEnd);
        
            % plot graph
            nodeLabels = arrayfun(@(n,i) sprintf('%d: %s',i,n.toStr()),obj.nodes,1:length(obj.nodes),'UniformOutput',false);
            han = plot(G,'Layout','layered','NodeLabel',nodeLabels);
            colormap([0 0 1]);
            drawnow;
        end
        
        % return handle only if desired
        if nargout == 0
            clear han;
        end
    end
end

methods (Access = protected)
    % copy the object
    function cp = copyElement(obj)
        cp = copyElement@matlab.mixin.Copyable(obj);
        for n = 1:length(obj.nodes)
            cp.nodes(n) = copy(obj.nodes(n));
        end
    end
end

methods (Access = private)
    % propagate the observations through the graph
    function propagate(obj,makeValid)
        obj.propagated = false(1,length(obj.nodes));
        obj.propagateRec(makeValid,1);
    end

    % recursive propagation helper function
    function propagateRec(obj,makeValid,idx)
        % first propagate the children (if they were not yet propagated) ...
        node = obj.nodes(idx);
        for c = setdiff(node.children',find(obj.propagated))
            obj.propagateRec(makeValid,c);
        end
        % ... and now the node at idx
        node.propagate(obj.nodes(node.children),makeValid);
        obj.propagated(idx) = true;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
