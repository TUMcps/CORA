classdef onlineAnalyzerNode < matlab.mixin.Copyable
% onlineAnalyzerNode - internal data structure of onlineReachSetAnalyzer
%                      represents a subformula of a given STL formula
%
% Syntax:
%    node = onlineAnalyzerNode(type,idx,parents)
%    node = onlineAnalyzerNode('AP',idx,parents,ap)
%    node = onlineAnalyzerNode('until',idx,parents,interval)
%
% Inputs:
%    type - the STL operator represented by the node
%          'AP' and 'until' are special types requiring additional arguments
%    idx - index in the overall graph structure of the onlineReachSetAnalyzer
%    parents - indices of the parents in the graph
%    ap - geometric representation of the associated atomic proposition
%    interval - time interval of the until operator
%
% Outputs:
%    node - generated onlineAnalyzerNode object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: onlineReachSetAnalyzer

% Authors:       Florian Lercher
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % basic properties
    type string
    signal fourValuedSignal
    prop (:,1)
    dropPropagated = true;
    
    % graph structure
    idx {mustBeInteger,mustBeNonnegative,mustBeScalarOrEmpty} = [];
    parents (:,1) {mustBeInteger,mustBeNonnegative} = [];
    children (:,1) {mustBeInteger,mustBePositive} = [];
    
    % properties only for until nodes
    interval {mustBeScalarOrEmpty} = [];
    
    % properties only for AP nodes
    ap {mustBeScalarOrEmpty} = [];
    tentative tentativeKleeneSignal = tentativeKleeneSignal.empty;
end

methods
    % constructor
    function obj = onlineAnalyzerNode(type,idx,parents,dropPropagated,varargin)

        inputArgsCheck({{type,'str',{'true','false','AP','~','|','&','until'}}});
        
        obj.type = type;
        obj.idx = idx;
        obj.parents = parents;
        obj.children = [];
        obj.prop = repmat(stlInterval(),size(parents));
        obj.dropPropagated = dropPropagated;
        obj.signal = aux_initSignal(type);

        % process additional arguments
        numBasicArgs = 4;
        switch type
            case {'true','false','~','|','&'}
                narginchk(0,numBasicArgs);
            case 'AP'
                narginchk(numBasicArgs+1,numBasicArgs+1);
                obj.ap = varargin{1};
                obj.tentative = tentativeKleeneSignal.emptySignal();
            case 'until'
                narginchk(numBasicArgs+1,numBasicArgs+1);
                obj.interval = varargin{1};
            otherwise
                assert(false,'Unreachable');
        end
    end

    % ------graph structure------
    % register a new parent in the graph structure
    function addParent(obj,parent)
        obj.parents = [obj.parents;parent];
        obj.prop = [obj.prop;stlInterval()];
    end

    % register a new child in the graph structure
    function addChild(obj,child)
        obj.children = [obj.children;child];
    end

    % ------observation------
    % observe a new part of the reachable set and update satisfaction signals
    function observeSet(obj,set,interval,loc,tentative)
        if ~ismember(obj.type,{'true','false','AP'})
            throw(CORAerror('CORA:notSupported','Only leaf nodes support observeSet'));
        end
        [lb,lc] = infimum(interval);
        for i = 1:length(obj.prop)
            [pUb,pRc] = supremum(obj.prop(i));
            if ~isempty(pUb) && (pUb > lb || (withinTol(pUb,lb,eps) && pRc && lc))
                throw(CORAerror('CORA:specialError','New observation for already propagated time is not allowed'));
            end
        end

        switch obj.type
            case {'true','false'}
                % do nothing, as the signal is already properly initialized
            case 'AP'
                canBeTrue = obj.ap.canBeTrue(set,loc);
                canBeFalse = obj.ap.canBeFalse(set,loc);
                
                if canBeTrue && ~canBeFalse
                    val = kleene.True;
                elseif ~canBeTrue && canBeFalse
                    val = kleene.False;
                elseif canBeTrue && canBeFalse
                    val = kleene.Unknown;
                else % ~canBeTrue && ~canBeFalse
                    throw(CORAerror('CORA:notDefined','AP can neither be true nor false'));
                end
                if tentative
                    % store the observation for later
                    obj.tentative = obj.tentative.observe(interval,val);
                else
                    % immediately update the satisfaction signal
                    obj.signal = obj.signal.set(interval,fourValued.fromKleene(val));
                end
            otherwise
                assert(false,'Unreachable');
        end
    end

    % ------propagation------
    % propagate observations from children
    function propagate(obj,childNodes,makeValid)
        switch obj.type
            case {'true','false'}
                % do nothing, as we have no children and no tentative observations
            case 'AP'
                % update the satisfaction signal with the tentative
                % observations during the interval makeValid
                if ~isemptyobject(makeValid)
                    final = obj.tentative.getInterval(makeValid);
                    obj.signal = obj.signal.merge(final);
                    obj.tentative = obj.tentative.setInvalid(makeValid);
                end
                % do nothing, as we have no children
            otherwise
                % compute the updated satisfaction signal from the children
                delta = obj.combine(childNodes);

                % merge the new signal with the old one
                obj.signal = obj.signal.merge(delta);

                % determine the conclusive interval
                % ignoring the times that have been propagated to all parents
                concInt = conclusiveInterval(obj.signal,obj.getAllProp());

                % set the propagated interval for all children
                % and drop the irrelevant prefixes of their signals
                for i = 1:length(childNodes)
                    childNodes(i).setPropagated(obj.idx,concInt);
                end
        end
    end

    % update the already propagated time interval
    function setPropagated(obj,parent,int)
        i = obj.parents == parent;
        if ~(obj.prop(i) == int) && contains(obj.prop(i),int)
            throw(CORAerror('CORA:specialError','New propagated is not allowed to be smaller than the old one'));
        end
        obj.prop(i) = int;

        if obj.dropPropagated % this should only be turned off for visualization purposes
            % drop the part of the signal that has been propagated to all parents
            obj.signal = obj.signal.set(obj.getAllProp(),fourValued.Inconclusive);
        end
    end

    % ------plotting------
    % convert to string for plots
    function str = toStr(obj)
        str = obj.type;
        if strcmp(obj.type,'until')
            str = sprintf('%s %s',str,obj.interval.toStr());
        end
    end
end

methods (Access = private)
    % combine signals of children according to the type
    function sig = combine(obj,childNodes)
        switch obj.type
            case {'true','false','AP'}
                throw(CORAerror('CORA:notSupported','Only non-leaf nodes support combine'));
            case '~'
                assert(length(childNodes) == 1);
                sig = ~childNodes(1).signal;
            case '&'
                assert(~isempty(childNodes));
                % childNodes(1).signal only determines the correct class here
                % it is not automatically passed to the static and_
                sig = childNodes(1).signal.and_(childNodes.signal);
            case '|'
                assert(~isempty(childNodes));
                % childNodes(1).signal only determines the correct class here
                % it is not automatically passed to the static or_
                sig = childNodes(1).signal.or_(childNodes.signal);
            case 'until'
                assert(length(childNodes) == 2);
                sig = childNodes(1).signal.until(obj.interval,childNodes(2).signal);
            otherwise
                assert(false,'Unreachable');
        end
    end

    % find the time interval that is already propagated to all parents
    function allProp = getAllProp(obj)
        ub = inf;
        rc = true;
        for i=1:length(obj.prop)
            if isemptyobject(obj.prop(i))
                % if we found a parent that did not receive any propagations, we can stop and return the empty interval
                ub = 0;
                rc = false;
                break
            else
                % otherwise, check if the propagation point is smaller than the current one
                [nUb, nRc] = supremum(obj.prop(i));
                if nUb < ub || withinTol(nUb,ub,eps) && ~nRc
                    ub = nUb;
                    rc = nRc;
                end
            end
        end
        allProp = stlInterval(0,ub,true,rc);
    end
end
end


% Auxiliary functions -----------------------------------------------------

function sig = aux_initSignal(type)
switch type
    case 'true'
        sig = fourValuedSignal.uniformSignal(fourValued.True);
    case 'false'
        sig = fourValuedSignal.uniformSignal(fourValued.False);
    otherwise
        sig = fourValuedSignal.uniformSignal(fourValued.Inconclusive);
end
end

% ------------------------------ END OF CODE ------------------------------
