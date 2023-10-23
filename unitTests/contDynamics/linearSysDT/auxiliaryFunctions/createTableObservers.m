function createTableObservers
% createTableObservers - creates table for the performance of set-based observers
%
% Syntax:
%    createTableObservers()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Matthias Althoff
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'linearSysDT' filesep 'results'];

% load file
name = 'TankLin_30_States_ZOrder_100_2021March29_10_59_24';
load(name);

% set format
formatSpec = '%#4.4g';

%% absolute values
% create file for table
fid = fopen([path filesep 'table_' name '_absolute.m'],'w');

% loop through observers
for iObserver = 1:length(estSet)
    % currentObserver
    currObserver = estSet{iObserver};
    % write name
    fprintf(fid, '%-13s', [currObserver.Name,' & ']);
    % write type
    fprintf(fid, '%-17s', [aux_determineType(currObserver),' & ']);
    % add ready for control option
    if iObserver <= 10
        fprintf(fid, '%-15s', '\xmark & ');
    else
        fprintf(fid, '%-15s', '\cmark & ');
    end
    % add radii
    radVec = currObserver.Performance.IRadius;
    for iRad = 1:length(radVec)
        fprintf(fid, '%s', [num2str(radVec(iRad), formatSpec),' & ']);
    end
    % add average
    fprintf(fid, '%s', [num2str(mean(radVec), formatSpec),' & ']);
    % print computation time in ms
    fprintf(fid, '%s\n', [num2str(1e3*currObserver.tIteration, formatSpec),' \\']);
    % collect values for relative table
    radMatrix(iObserver,:) = radVec;
end

%close file
fclose(fid);


%% relative values
% create file for table
fid = fopen([path filesep 'table_' name '_relative.m'],'w');

% find best radius for each dimension
bestRad = min(radMatrix);

% rompute relative values for radMatrix
for iCol = 1:length(radMatrix(1,:))
    radMatrixRel(:,iCol) = radMatrix(:,iCol)/bestRad(iCol);
end

% loop through observers
for iObserver = 1:length(estSet)
    % currentObserver
    currObserver = estSet{iObserver};
    % write name
    fprintf(fid, '%-13s', [currObserver.Name,' & ']);
    % write type
    fprintf(fid, '%-17s', [aux_determineType(currObserver),' & ']);
    % add ready for control option
    if iObserver <= 10
        fprintf(fid, '%-15s', '\xmark & ');
    else
        fprintf(fid, '%-15s', '\cmark & ');
    end
    % add radii
    radVecRel = radMatrixRel(iObserver,:);
    for iRad = 1:length(radVecRel)
        fprintf(fid, '%s', [num2str(radVecRel(iRad), formatSpec),' & ']);
    end
    % add average
    fprintf(fid, '%s', [num2str(mean(radVecRel), formatSpec),' & ']);
    % print computation time in ms
    fprintf(fid, '%s\n', [num2str(1e3*currObserver.tIteration, formatSpec),' \\']);
    
    %% add special rows at the end of each observer type
    if iObserver == 10
        fprintf(fid, '%s\n', '\midrule \multicolumn{11}{c}{set-propagation observers} \\ \midrule');
    elseif iObserver == 15
        fprintf(fid, '%s\n', '\midrule \multicolumn{11}{c}{interval observer} \\ \midrule');
    end
end
%% last row
% print seperator
fprintf(fid, '%s\n', '\midrule \multicolumn{11}{c}{smallest absolute radii} \\ \midrule');
% new columns
fprintf(fid, '%-13s', ' & & & ');
% best values
for iRad = 1:length(bestRad)
    fprintf(fid, '%s', [num2str(bestRad(iRad), formatSpec),' & ']);
end


%close file
fclose(fid);


%% relative values without average
% create file for table
fid = fopen([path filesep 'table_' name '_relativeNoAverage.m'],'w');

% find best radius for each dimension
bestRad = min(radMatrix);

% rompute relative values for radMatrix
for iCol = 1:length(radMatrix(1,:))
    radMatrixRel(:,iCol) = radMatrix(:,iCol)/bestRad(iCol);
end

% loop through observers
for iObserver = 1:length(estSet)
    % currentObserver
    currObserver = estSet{iObserver};
    % write name
    fprintf(fid, '%-13s', [currObserver.Name,' & ']);
    % write type
    fprintf(fid, '%-17s', [aux_determineType(currObserver),' & ']);
    % add ready for control option
    if iObserver <= 10
        fprintf(fid, '%-15s', '\xmark & ');
    else
        fprintf(fid, '%-15s', '\cmark & ');
    end
    % add radii
    radVecRel = radMatrixRel(iObserver,:);
    for iRad = 1:length(radVecRel)
        fprintf(fid, '%s', [num2str(radVecRel(iRad), formatSpec),' & ']);
    end
    % print computation time in ms
    fprintf(fid, '%s\n', [num2str(1e3*currObserver.tIteration, formatSpec),' \\']);
    
    %% add special rows at the end of each observer type
    if iObserver == 10
        fprintf(fid, '%s\n', '\midrule \multicolumn{10}{c}{set-propagation observers} \\ \midrule');
    elseif iObserver == 15
        fprintf(fid, '%s\n', '\midrule \multicolumn{10}{c}{interval observer} \\ \midrule');
    end
end
%% last row
% print seperator
fprintf(fid, '%s\n', '\midrule \multicolumn{10}{c}{smallest absolute radii} \\ \midrule');
% new columns
fprintf(fid, '%-13s', ' & & & ');
% best values
for iRad = 1:length(bestRad)
    fprintf(fid, '%s', [num2str(bestRad(iRad), formatSpec),' & ']);
end


%close file
fclose(fid);


end


% Auxiliary functions -----------------------------------------------------

% determine type of observer
function type = aux_determineType(observer)
    % obtain first set
    firstSet = observer.EstStates.timePoint.set{1};
    % zonotope?
    if isa(firstSet,'zonotope')
        type = 'zonotope';
    elseif isa(firstSet,'ellipsoid')
        type = 'ellipsoid';
    elseif isa(firstSet,'conZonotope')
        type = 'constr. zono.';
    end
        
end

% ------------------------------ END OF CODE ------------------------------
