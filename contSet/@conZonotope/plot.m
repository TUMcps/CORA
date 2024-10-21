function han = plot(cZ,varargin)
% plot - plots a projection of a constrained zonotope
%
% Syntax:
%    han = plot(cZ)
%    han = plot(cZ,dims)
%    han = plot(cZ,dims,type)
%
% Inputs:
%    cZ - conZonotope object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'Splits',splits> - number of splits for refinement
%               <'Template',dirs> - template directions
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    figure;
%    plot(cZ,[1,2],'r');
%
%    figure;
%    plot(cZ,[1,2],'FaceColor','b','Splits',4);
%
%    figure;
%    plot(cZ,[1,2],'r','Template',16);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       11-May-2018
% Last update:   15-July-2020 (MW, merge with plotFilled|Template|Split)
%                25-May-2022 (TL, 1D Plotting)
%                16-December-2022 (MW, add iterative method for 2D plots)
%                05-April-2023 (TL, clean up using plotPolygon)
%                27-April-2023 (VG, check if cZ is feasible)
%                09-May-2023 (TL, bugfix split plotting)
%                15-October-2024 (TL, use contSet/plot for default mode)
% Last revision: 12-July-2023 (TL, restructure)
%                15-October-2024 (TL, split into plot1D/plot2D/plot3D for default plotting)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[cZ,dims,NVpairs,mode,splits,numDir] = aux_parseInput(cZ,varargin{:});

% 2. preprocess
[cZ,dims] = aux_preprocess(cZ,dims);

% 3. plot n-dimensional set
if mode == 1
    % default plotting
    han = plot@contSet(cZ,dims,NVpairs{:});
else
    % special plotting
    han = aux_plotNd(cZ,dims,NVpairs,mode,splits,numDir);
end

% 4. clear han
if nargout == 0
    clear han
end

end


% Auxiliary functions -----------------------------------------------------

function [cZ,dims,NVpairs,mode,splits,numDir] = aux_parseInput(cZ,varargin)
    % parse input

    % default settings
    dims = setDefaultValues({[1,2]},varargin);

    % check input args
    inputArgsCheck({{cZ, 'att', 'conZonotope'},
        {dims,'att','numeric',{'nonempty','vector','integer','positive'}}})

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    
    % process linespec and Name-Value pairs
    NVpairs = readPlotOptions(varargin(2:end));
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
    [NVpairs,numDir] = readNameValuePair(NVpairs,'Template','isscalar');
    
    % choose plot mode
    if ~isempty(splits) && ~isempty(numDir)
        % error if both given
        throw(CORAerror('CORA:specialError','Choose either Splits or Template.'));
    elseif ~isempty(splits)
        % plot mode 'split'
        mode = 2;
    elseif ~isempty(numDir)
        % plot mode 'template'
        mode = 3;
    else
        % default plot mode 'standard'
        mode = 1;
    end
end

function [cZ,dims] = aux_preprocess(cZ,dims)
    % project the object to the N=dimensional subspace
    cZ = project(cZ,dims);
    dims = 1:length(dims);
end

function han = aux_plotNd(cZ,dims,NVpairs,mode,splits,numDir)
    % check if constraints are feasible
    if representsa_(cZ,'emptySet',eps)
        % plot empty set
        han = plot(emptySet(numel(dims)), dims, NVpairs{:});
        return
    end

    % plot modes: standard (1), template (2), splits (3)
    switch mode
        case 2
            han = aux_plotSplit(cZ,splits,dims,NVpairs);
        case 3
            V = vertices(cZ,'template',numDir);
            han = plotPolygon(V,'ConvHull',true,NVpairs{:});
        otherwise
            % default plot mode
            han = plot@contSet(cZ,dims,NVpairs{:});
    end

end

function han = aux_plotSplit(cZ,splits,dims,NVpairs)

    % recursively split the constrained zonotope
    cZSplit = {cZ};

    for i = 1:splits
        % init
        listTemp = cell(length(cZSplit)*2,1);
        counter = 1;

        % loop over all sets at the current recursion level
        for j = 1:length(cZSplit)

           % calculate radius of the interval over-approximation as a
           % heuristic indicating which dimension should be best splitted
           inter = interval(cZSplit{j});
           r = rad(inter);

           % split the set
           [~,ind] = max(r);
           temp = split(cZSplit{j},ind);

           % update variables
           listTemp{counter} = temp{1};
           listTemp{counter+1} = temp{2};
           counter = counter + 2;
        end

        cZSplit = listTemp;
    end

    % convert splitted sets to intervals
    Is = cell(1,length(cZSplit));
    for i=1:length(cZSplit)
        Is{i} = zonotope(cZSplit{i});
    end

    % plot all sets as one
    han = plotMultipleSetsAsOne(Is,dims,NVpairs);
end

% ------------------------------ END OF CODE ------------------------------
