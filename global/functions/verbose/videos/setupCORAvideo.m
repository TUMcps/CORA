function [vidObj,fig] = setupCORAvideo(filename,varargin)
% setupCORAvideo - prepares everything to record a CORA video
%
% Syntax:
%    res = setupCORAvideo(filename,varargin)
%
% Inputs:
%    filename - char, video file name (with extension *.mp4)
%    varargin - name-value-pairs
%       'Title' - title of video
%       'Description' - description (also controls layout)
%       'XLabel' - xlabel
%       'XLim' - xlim
%       'YLabel' - ylabel
%       'YLim' - ylim
%       'LegendLocation' - location of legend
%       'FrameRate' - frame rate of video
%
% Outputs:
%    vidObj - VideoWriter
%    fig - figure
%
% See also:
%    CORAvideo_snippets

% Authors:       Tobias Ladner
% Written:       05-April-2024
% Last update:   11-December-2024 (TL, split setup and recording)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp('Setting up CORA video:')

% parse input ---

[XLabel,XLim,YLabel,YLim,ZLabel,ZLim,Title,Description,LegendLocation,FrameRate] = aux_parseInput(varargin);

% setup ---

% set up figure
fig = aux_setUpFigure(XLabel,XLim,YLabel,YLim,ZLabel,ZLim,Title,Description,LegendLocation);

% set up video
vidObj = VideoWriter(filename, 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = FrameRate;
open(vidObj);

% display path
videoName = [vidObj.Path filesep vidObj.Filename];
fprintf('- Path: <a href="matlab:winopen(''%s'');">%s</a> \n', vidObj.Path, videoName)

end


% Auxiliary functions -----------------------------------------------------

function [XLabel,XLim,YLabel,YLim,ZLabel,ZLim,Title,Description,LegendLocation,FrameRate] = aux_parseInput(NVpairs)
    
    % read input args
    [NVpairs,XLabel] = readNameValuePair(NVpairs,'XLabel',{},'x_{(1)}');
    [NVpairs,XLim] = readNameValuePair(NVpairs,'XLim');
    [NVpairs,YLabel] = readNameValuePair(NVpairs,'YLabel',{},'x_{(2)}');
    [NVpairs,YLim] = readNameValuePair(NVpairs,'YLim');
    [NVpairs,ZLabel] = readNameValuePair(NVpairs,'ZLabel',{},'x_{(3)}');
    [NVpairs,ZLim] = readNameValuePair(NVpairs,'ZLim');
    [NVpairs,Title] = readNameValuePair(NVpairs,'Title',{},'');
    [NVpairs,Description] = readNameValuePair(NVpairs,'Description',{},'');
    [NVpairs,LegendLocation] = readNameValuePair(NVpairs,'LegendLocation',{},'northeast');
    [NVpairs,FrameRate] = readNameValuePair(NVpairs,'FrameRate',{},30);
 
end

function fig = aux_setUpFigure(XLabel,XLim,YLabel,YLim,ZLabel,ZLim,Title,Description,LegendLocation)
    
    % set up figure
    fig = figure; 
    set(fig,'WindowStyle','normal');
    width = 2160;
    height = 2160;

    % background color white
    backgroundColor = [1 1 1];
    set(gcf, "Color", backgroundColor);

    % font size
    fontSize = 36;
    fontSizeTitle = ceil(fontSize * 1.33);
    fontSizeSmall = ceil(fontSize / 1.33);

    if isempty(Description)
        % just show plot
        aux_setFigureSize(fig, width, height)

        % and write title above plot
        hold on; box on;
        titleHandle = title(Title);
        titleHandle.FontSize = fontSizeTitle;
    else
        % make space for title and description
        width = 3840;
        aux_setFigureSize(fig, width, height)

        % plot title
        dim = [.025 .8 .45 .1]; 
        annotation('textbox',dim,'String',Title,'EdgeColor','none','FontWeight','bold','FontSize',fontSizeTitle)
        
        % plot description
        dim = [.025 .3 .45 .5]; 
        annotation('textbox',dim,'String',compose(Description),'EdgeColor','none','FontSize',fontSize)

        % show image
        img = imread('coraLogo.png','BackgroundColor',[1 1 1]);
        s = size(img);
        imgwidth = 0.075;
        axCORA = axes( ...
            'Units','pixels', ...
            'OuterPosition',[width*0.025 height*0.05 width*imgwidth s(1)/s(2)*imgwidth*width+fontSizeSmall/0.75] ...
        );
        image(axCORA,img)
        axCORA.Toolbar.Visible = 'off';
        disableDefaultInteractivity(axCORA)
        set(axCORA,'xtick',[])
        set(axCORA,'ytick',[])
        set(axCORA,'XColor',backgroundColor,'YColor',backgroundColor,'TickDir','out')
        website = xlabel(axCORA, 'https://cora.in.tum.de');
        website.Color = 'k';
        website.FontSize = fontSizeSmall;

        % create subplot on the right
        subplot(1,2,2); hold on; box on;
    end

    fprintf('- Quality: %ix%i\n', width, height)

    % set up axis
    ax = gca();
    disableDefaultInteractivity(ax)
    ax.Toolbar.Visible = 'off';
    set(ax, "FontSize", fontSize)
    
    % label and limits
    if ~isempty(XLabel)
        xlabel(XLabel); 
    end
    if ~isempty(XLim)
        xlim(XLim)
    end
    if ~isempty(YLabel)
        ylabel(YLabel); 
    end
    if ~isempty(YLim)
        ylim(YLim)
    end
    if ~isempty(ZLabel)
        zlabel(ZLabel); 
    end
    if ~isempty(ZLim)
        zlim(zLim)
    end

    % set up legend
    legend('Location',LegendLocation)
end

function aux_setFigureSize(fig, width, height)
    screenSize = get(0, 'ScreenSize');
    left = (screenSize(3) - width) / 2;
    bottom = 0;
    set(fig, 'OuterPosition', [left, bottom, width, height]);
end

% ------------------------------ END OF CODE ------------------------------
