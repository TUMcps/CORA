function han = plot(R,varargin)
% plot - plots a projection of the reachable set
%
% Syntax:
%    han = plot(R)
%    han = plot(R,dims)
%    han = plot(R,dims,plotOptions)
%
% Inputs:
%    R - reachSet object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs),
%        including added pairs:
%          'Order', <order> (polynomial) zonotope order
%          'Splits', <splits> number of splits for refinement
%          'Unify', <true/false> compute union of all reachable sets
%          'Set', <whichset> corresponding to
%                   ti ... time-interval reachable set (default)
%                   tp ... time-point reachable set (default if no ti)
%                   y  ... time-interval algebraic set
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-June-2020
% Last update:   15-July-2020 (MW, merge with plotFilled, plotOptions)
%                29-October-2021 (MW, remove linespec, 'Filled')
%                11-July-2023 (VG, bug fix unify last set not plotted)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------


% 1. parse input
[R,dims,NVpairs,order,splits,unify,totalsets,whichset] = aux_preprocess(R,varargin{:});

% 2. plot reachable sets
han = aux_plotReachSet(R,dims,NVpairs,order,splits,unify,totalsets,whichset);

% 3. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [R,dims,NVpairs,order,splits,unify,totalsets,whichset] = aux_preprocess(R,varargin)
    % parse input

    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{R,'att','reachSet','nonempty'};
                    {dims,'att','numeric',{'nonempty','vector','integer','positive'}}});
    
    % check name-value pairs
    NVpairs = readPlotOptions(varargin(2:end),'reachSet');
    [NVpairs,order] = readNameValuePair(NVpairs,'Order','isscalar');
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
    [NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',false);
    [NVpairs,totalsets] = readNameValuePair(NVpairs,'UnifyTotalSets','isscalar',1);
    [NVpairs,whichset] = readNameValuePair(NVpairs,'Set','ischar','ti');
    
    % check dimension
    if length(dims) < 2
        throw(CORAerror('CORA:plotProperties',2));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    elseif length(dims) == 3
        unify = false;
    end
    
    % check which set has to be plotted
    whichset = aux_checkSet(R,whichset);
end

function res = aux_hasTimeInterval(R)
    % check whether reachSet has a time interval set
    res = ~isempty(R(1).timeInterval);
end

function whichset = aux_checkSet(R,whichset)

% must be character vector for switch-expression to work properly
if isempty(whichset)
    % default value
    if ~isempty(R(1).timeInterval)
        whichset = 'ti';
    else
        whichset = 'tp';
    end
end

switch whichset
    case 'ti'
        if ~aux_hasTimeInterval(R) && ~isempty(R(1).timePoint)
            warning("No time-interval reachable set. Time-point reachable set plotted instead.");
            whichset = 'tp';
        end
        
    case 'tp'
        % no issues (should always be there)

    case 'y'
        if isempty(R(1).timeInterval.algebraic)
            throw(CORAerror('CORA:emptyProperty'));
        end

    otherwise
        % default value
        if isempty(whichset)
            whichset = 'ti';
            if isempty(R(1).timeInterval) 
                if ~isempty(R(1).timePoint)
                    whichset = 'tp';
                else
                    throw(CORAerror('CORA:emptySet'));
                end
            end
        else
            % has to be one of the above three cases...
            throw(CORAerror('CORA:wrongValue','Set',...
                'has to be ''ti'', ''tp'' or ''y''.'));
        end
end

end

function han = aux_plotReachSet(R,dims,NVpairs,order,splits,unify,totalsets,whichset)
    % plots the reachable set
    
    % check if the reachable sets should be unified to reduce the storage size
    % of the resulting figure
    if unify
        % unify and plot the reachable sets
        han = aux_plotUnified(R,dims,NVpairs,order,splits,totalsets,whichset);
      
    else
        % plot reachable sets individually
        han = aux_plotSingle(R,dims,NVpairs,order,splits,whichset);
    end

end

function han = aux_plotUnified(R,dims,NVpairs,order,splits,totalsets,whichset)
    % unify and plot the reachable sets

    % number of all sets that are to be unified
    nrAllSets = aux_nrAllSets(R,whichset);
    idxSplit = splitIntoNParts(nrAllSets,totalsets);
    idxCurr = 0;

    % init
    pgons = {};
    pgon = polygon();
    warOrig = warning;
    warning('off','all');
    
    % loop over all reachable sets
    for i = 1:size(R,1)
        
        % get desired set
        switch whichset
            case 'ti'
                Rset = R(i,1).timeInterval.set;
            case 'tp'
                Rset = R(i,1).timePoint.set;
            case 'y'
                Rset = R(i,1).timeInterval.algebraic;
        end
        
        % loop over all time steps
        for j = 1:length(Rset)

            % increment index
            idxCurr = idxCurr + 1;

            % project set to desired dimensions
            temp = project(Rset{j},dims);

            % order reduction
            if ~isempty(order)
                temp = reduce(temp,'girard',order);
            end
           
            % convert set to polygon object
            if isa(temp,'polyZonotope') || isa(temp,'conPolyZono')
                if isempty(splits)
                    V = vertices(zonotope(temp));
                    temp = polygon(V(1,:),V(2,:));
                else
                    temp = polygon(temp,splits);   
                end
            elseif isa(temp,'zonotope') || isa(temp,'interval') || ...
                   isa(temp,'polytope') || isa(temp,'conZonotope')
                V = vertices(temp);
                temp = polygon(V(1,:),V(2,:));
            elseif isa(temp,'zonoBundle')                
                % compute polytopes
                P = cell(temp.parallelSets,1);
                for p = 1:temp.parallelSets
                    % delete zero-length generators
                    Z = compact_(temp.Z{p},'zeros',eps);
                    % convert to polyshape (Matlab built-in class)
                    temp_poly = polygon(Z);
                    V = temp_poly(:,2:end);
                    P{p} = polyshape(V(1,:),V(2,:));
                end

                % intersect polytopes
                Pint = P{1};
                for p = 2:temp.parallelSets
                    Pint = intersect(Pint,P{p});
                end

                % get vertices, init polygon
                V = Pint.Vertices';
                temp = polygon(V(1,:),V(2,:));
            else
                warning('CORA: Setting "Unify" is not supported for this set representation (%s)! Plotting them individually instead.', class(temp));
                han = aux_plotSingle(R,dims,NVpairs,order,splits,whichset);
                return
            end
            
            % unite all polygons
            pgon = pgon | temp;

            if any(idxCurr == idxSplit(:,2))
                % end of partition reached -> plot
                pgons{end+1} = pgon;
                
                % reset pgon
                pgon = polygon();
            end
        end
    end
    
    warning(warOrig);

    % plot polygons
    han = plotMultipleSetsAsOne(pgons,[1,2],NVpairs);

end

function nrAllSets = aux_nrAllSets(R,whichset)
    % compute number of all sets that are to be plotted (only: 'Unify',true)
    % we will then partition this number into the desired number of sets given
    % by the name-value pair 'UnifyTotalSets',<totalsets> and plot <totalsets>
    % different sets (to avoid increasingly time-consuming '|'-operation of
    % polygons)
    
    % init total number
    nrAllSets = 0;
    
    % all branches of the reachSet object
    for i=1:size(R,1)
        % choose correct set
        switch whichset
            case 'ti'
                addSets = size(R(i,1).timeInterval.set,1);
            case 'tp'
                addSets = size(R(i,1).timePoint.set,1);
            case 'y'
                addSets = size(R(i,1).timeInterval.algebraic,1);
        end
        nrAllSets = nrAllSets + addSets;
    end
end

function han = aux_plotSingle(R,dims,NVpairs_base,order,splits,whichset)
    % plot reachable sets individually

    % init
    sets = {};
    NVpairs = {};

    % iterate over all reachable sets
    for i = 1:size(R,1)

        % get desired set
        switch whichset
            case 'ti'
                Rset = R(i,1).timeInterval.set;
            case 'tp'
                Rset = R(i,1).timePoint.set;
            case 'y'
                Rset = R(i,1).timeInterval.algebraic;
        end

        % loop over all time steps
        for j = 1:length(Rset)

            % project set to desired dimensions
            R_proj = project(Rset{j},dims);

            % order reduction
            if ~isempty(order)
                R_proj = reduce(R_proj,'girard',order);
            end

            % plot the set
            if isa(R_proj,'polyZonotope')
                if isempty(splits)
                    sets{end+1} = zonotope(R_proj);
                    NVpairs{end+1} = NVpairs_base;
                else
                    sets{end+1} = R_proj;
                    NVpairs{end+1} = [{'Splits',splits}, NVpairs_base];
                end
            else
                sets{end+1} = R_proj;
                NVpairs{end+1} = NVpairs_base;
            end
        end
    end

    % plot sets
    han = plotMultipleSetsAsOne(sets,1:length(dims),NVpairs);

end

% ------------------------------ END OF CODE ------------------------------
