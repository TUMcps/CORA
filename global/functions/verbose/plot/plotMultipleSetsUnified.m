function han = plotMultipleSetsUnified(sets,varargin)
% plotMultipleSetsUnified - plots a list of sets to the current figure
%    while attempting to unify the resulting polygon objects
%    (this results in faster plotting and helps when exporting to tikz)
%
% Syntax:
%    han = plotMultipleSetsUnified(sets,dims,NVpairs)
%    han = plotMultipleSetsUnified(sets,dims,NVpairs,purpose)
%
% Inputs:
%    sets - cell array of sets
%    dims - desired dimensions to plot
%    NVpairs - plot settings name-value pairs
%        'UnifyTotalSets', <numTotalSets> - number of total unified sets
%        'Splits', <splits> - number of splits for refinement
%        'Order', <order> - (polynomial) zonotope order
%    purpose - purpose - information about plot
%           or cell array containing purpose for each set
%
% Outputs:
%    han - handle to the graphics objects
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet/plot, plotMultipleSetsAsOne

% Authors:       Tobias Ladner
% Written:       11-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[dims,NVpairs,purpose,numTotalSets,idxSplit,order,splits] = aux_parseInput(sets,varargin{:});

% init
pgons = cell(1,min([numTotalSets,numel(sets)]));
pgon = polygon();
warOrig = warning;
warning('off','all');

% unify sets
for i=1:numel(sets)

    % read out set
    S_i = sets{i};

    % project
    S_i = project(S_i,dims);

    % order reduction
    if ~isempty(order)
        S_i = reduce(S_i,'girard',order);
    end
    
    % convert set to polygon object
    try
        % special treatment for polynomial zonotopes
        if isa(S_i,'polyZonotope') || isa(S_i,'conPolyZono')
            % consider number of splits
            if isempty(splits)
                pgon_i = polygon(zonotope(S_i));
            else
                pgon_i = polygon(S_i,splits);   
            end
        else
            % convert directly to polygon
            pgon_i = polygon(S_i);
        end
    catch ME
        CORAwarning('CORA:plot','Setting "Unify" failed for %s (Error: %s)!\nPlotting them individually instead.', class(S_ij),ME.message);
        han = plotMultipleSetsAsOne(sets,dims,NVpairs,purpose);
        return
    end

    % unify
    timerVal = tic;
    pgon = pgon | pgon_i;
    duration = toc(timerVal);

    % choosing small threshold as it quickly adds up with many plotted sets
    durationThreshold = 0.001; % seconds

    % reset unification either if on ...
    if i == numel(sets) ... % (i) last iteration
        ... % (ii) timer-based
        || (isempty(numTotalSets) && duration > durationThreshold) ... 
        ... % (iii) pre-determined splits
        || (~isempty(numTotalSets) && any(i == idxSplit(:,2))) ... 
        % reset
        pgons{i} = pgon;
        pgon = polygon();
    end
end

% filter empty objects
pgons = pgons(cellfun(@(pgon) isa(pgon,'polygon'),pgons,'UniformOutput',true));

% split up polygons with multiple regions
pgons = cellfun(@(pgon) pgon.getRegions()', pgons,'UniformOutput',false);
pgons = [pgons{:}];

% plot sets
han = plotMultipleSetsAsOne(pgons,1:2,NVpairs,purpose);

% clean up
warning(warOrig);
% clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [dims,NVpairs,purpose,numTotalSets,idxSplit,order,splits] = aux_parseInput(sets,varargin)
   % set default values
   [dims,NVpairs,purpose] = setDefaultValues({[1,2],{},'none'},varargin);
   % check input
   inputArgsCheck({ ...
       {sets,'att','cell'}
       {dims,'att','numeric',{'nonempty','vector','integer','positive'}}
       {NVpairs,'att','cell'}
   });

   % check name-value pairs
   [NVpairs,numTotalSets] = readNameValuePair(NVpairs,'UnifyTotalSets',{},[]); % [] means 'automatic'
   [NVpairs,order] = readNameValuePair(NVpairs,'Order','isscalar');
   [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');

   % determine indices which sets are unified
   idxSplit = splitIntoNParts(numel(sets),min([numel(sets),numTotalSets]));
end

% ------------------------------ END OF CODE ------------------------------
