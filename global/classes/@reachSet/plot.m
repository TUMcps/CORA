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
%          'Order', <order> - (polynomial) zonotope order
%          'Splits', <splits> - number of splits for refinement
%          'Unify', <true/false> - compute union of all reachable sets
%          'UnifyTotalSets', <numTotalSets> - number of total unified sets
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
%                10-April-2024 (TL, bug fix UnifyTotalSets too large)
%                17-December-2024 (TL, added option to plot empty objects)
%                11-July-2025 (TL, automatic unification)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------


% 1. parse input
[R,dims,NVpairs,unify,whichset] = aux_preprocess(R,varargin{:});

% 2. plot reachable sets
han = aux_plotReachSet(R,dims,NVpairs,unify,whichset);

% 3. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [R,dims,NVpairs,unify,whichset] = aux_preprocess(R,varargin)
    % parse input

    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{R,'att','reachSet','nonempty'};
                    {dims,'att','numeric',{'nonempty','vector','integer','positive'}}});
    
    % check name-value pairs
    NVpairs = readPlotOptions(varargin(2:end),'reachSet');
    [NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',numel(dims) == 2);
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

function whichset = aux_checkSet(R,whichset)

% must be character vector for switch-expression to work properly
if isempty(whichset)
    % default value
    if ~hasTimeInterval(R)
        whichset = 'ti';
    else
        whichset = 'tp';
    end
end

switch whichset
    case 'ti'
        if ~hasTimeInterval(R) && ~isempty(R(1).timePoint)
            CORAwarning('CORA:plot',"No time-interval reachable set. Time-point reachable set plotted instead.");
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

function han = aux_plotReachSet(R,dims,NVpairs,unify,whichset)
    % plots the reachable set

    % empty case
    if isemptyobject(R)
        han = plotPolygon(nan(numel(dims),0),NVpairs{:});
        return
    end

    % gather list of sets to plot
    sets = arrayfun(@(R) aux_gatherSetsScalar(R,whichset),R,'UniformOutput',false);
    sets = [sets{:}];

    % 3D plots and probZonotopes are handled differently
    if numel(dims) == 3 || any(cellfun(@(S) isa(S,'probZonotope'), sets))
        han = plotMultipleSetsAsOne(sets,dims,NVpairs);
    elseif unify
        % unify and plot the reachable sets
        han = plotMultipleSetsUnified(sets,dims,NVpairs);
    else
        % no unification:
        % use plotMultipleSetsUnified anyway due to pre-processing handled there
        % while hard-coding to not unify any sets
        han = plotMultipleSetsUnified(sets,dims,[NVpairs,{'UnifyTotalSets',numel(sets)}]);
    end

end

function sets = aux_gatherSetsScalar(R,whichset)
    % gathers all sets of a scalar reachable set based on given whichset    
    switch whichset
        case 'ti'
            sets = R.timeInterval.set;
        case 'tp'
            sets = R.timePoint.set;
        case 'y'
            sets = R.timeInterval.algebraic;
    end

    % make sure it's a row vector
    sets = reshape(sets,1,[]);
end

% ------------------------------ END OF CODE ------------------------------
