%% HOW TO RECORD VIDEOS IN CORA
%
% CORA enables recording Videos of everything you do directly in MATLAB.
% This files provides you with several code snippets to get you started.
% Generally, you should follow the following steps, which are also
% described in more detail in the individual code snippets below.
% _________________________________________________________________________
% Basic Steps to Record a CORA Video:
% -------------------------------------------------------------------------
% 0. Gather all required variables for the video
% -------------------------------------------------------------------------
% 1. Setup video and figure
%    -> setupCORAvideo(...)
% 2. Record your animations
%    a. Recording reachable set and simulations
%       -> recordReachableSet(...)
%    b. Animate sets as needed
%       -> animateFromTo(...)
% 3. Finish CORA video
%    -> finishCORAvideo(...)
% -------------------------------------------------------------------------
% 4. (Optional) Loop video and add music
%    see code snippet below
% 5. (Optional) Save title, description, and thumbnail
%    see code snippet below
% _________________________________________________________________________
%
% See also: 
%     setupCORAvideo, finishCORAvideo,
%     recordReachableSet, animateFromTo, 
%     writeFreezedFrames

% Authors:       Tobias Ladner
% Written:       12-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% 0. GATHER ALL REQUIRED VARIABLES FOR THE VIDEO -------------------------

% e.g. set breakpoint after everything is computed / at the end

%% 1. SETUP VIDEO AND FIGURE ----------------------------------------------

% settings ---

basepath = './CORAvideos/yy-mm-example-name';
exampleName = 'Discrete Time Systems';
mkdir(basepath)
filenameVideo = [basepath '/video.mp4'];

% title and description ---

% (grouped by topic, execute/adapt as needed)

%% Reachability Analysis for Continuous Systems

titleVideo = 'Reachability Analysis for Continuous Systems';

descVideo = "CORA computes reachable sets for linear systems, nonlinear systems, as well as for systems with constraints. " + ...
    "Continuous as well as discrete time models are supported. Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitly considered. " + ...
    "In addition, CORA also provides capabilities for the simulation of dynamical models.\n" + ...
    "\n" + ...
    "The example on the right formally verifies the stability of the Van-der-Pol oscillator system. " + ...
    "For details of the system, please visit Sec. 3.1 of the ARCH-COMP19 category report: Continuous and hybrid systems with nonlinear dynamics.";


%% Reachability Analysis for Hybrid Systems

titleVideo = 'Reachability Analysis for Hybrid Systems';

descVideo = "CORA is capable to calculate the reachable sets for hybrid systems. " + ...
    "All implemented dynamic system classes can be used to describe the different continuous flows for the discrete system states. " + ...
    "Further, multiple different methods for the calculation of the intersections with guard sets are implemented in CORA.\n" + ...
    "\n" + ...
    "The example on the right shows a ball bouncing off a table over time. " + ...
    "This benchmark is particularly difficult to verify due to the instant transitions every time the ball hits the table.";


%% Geometric Sets

titleVideo = 'Set-Based Computing';

descVideo = "CORA has a modular design, making it possible to use the capabilities of the various set representations for other purposes besides reachability analysis. " + ...
    "The toolbox implements vector set representation, e.g., intervals, zonotopes, polytopes, and Taylor models, as well as matrix set representations such as interval matrices and matrix zonotopes.\n" + ...
    "\n" + ...
    "The example on the right shows an example of a set-based algorithm. " + ...
    "The algorithm first propagates the initial set using affine maps and adds small error terms using the Minkowski sum for nine steps. " + ...
    "Finally, the union of two propagated sets and the intersection of the last two sets are computed.";


%% Formal Verification of Neural Networks

titleVideo = 'Formal Verification of Neural Networks';

descVideo = "CORA enables the formal verification of neural networks, both in open-loop as well as in closed loop scenarios. " + ...
    "Open-loop verification refers to the task where properties of the output set of a neural network are verified, e.g. correctly classified images given noisy input. " + ...
    "In closed-loop scenarios, the neural network is used as controller of a dynamic system and is neatly integrated in the reachability algorithms above, e.g. controlling a car while keeping a safe distance.\n" + ...
    "\n" + ...
    "The example on the right shows how the output of each layer of a small toy network with [2,3,2] neurons and sigmoid activation is enclosed. " + ...
    "Samples are shown with black dots, the enclosing set is shown in blue. " + ...
    "CORA ensures that all points from the input set are contained in the output enclosure of each layer and thus of the entire network.";


%% setup video and figure ---

[vidObj,fig] = setupCORAvideo(filenameVideo,'title',titleVideo,'Description',descVideo);

%% 2. RECORD YOUR ANIMATIONS ----------------------------------------------

%% 2.a Recording Reachable Set and Simulations (plot)

dims = 1:2;
recordReachableSet(vidObj,fig,'ReachSet',R,'SimResult',simRes,'Dimensions',dims, 'Specification', spec)

%% 2.a Recording Reachable Set and Simulations (plotOverTime)

dims = 1;
recordReachableSet(vidObj,fig, ...
    'ReachSet',R,'SimResult',simRes,'Dimensions',dims, ...
    'Title', titleVideo,'Description', descVideo)

%% 2.a Recording Reachable Set and Simulations (plot multiple reachable sets)

dims = 1;
recordReachableSet(vidObj,fig, ...
    'ReachSets',{ ...
        project(R_distances,1),'Distance',CORAcolor('CORA:safe'); ...
        project(R_distances,2),'Safe distance',CORAcolor('CORA:unsafe') ...
    },'SimResult',simResDistances,'Dimensions', dims);

%% 2.b Animate Sets as Needed

% Generally, the animation is defined via an animation struct specifying
% the start set/points (.fromS), the dimensions to plot (.dims), the
% name-value pairs for plotting (.fromNVpairs) and a transformation
% (.transformFun or .toS). 
% 
% Defining a transformation function @(S,t) allows you to specify how a
% set/points S is transformed over time t \in [0,1]. For examples, please 
% see the ./transformFuns folder. Alternatively, you can define a final
% set/points .toS and let CORA create a smooth transformation in between
% (convex combination over time via Minkowski sum; make sure that the 
% chosen set representation works well for that!).
%
% Finally, you might want to update the color and alpha values of your
% animated set/points over time. The (updated) name-value pairs can be
% specified in .toNVpairs. Here, you only need to specify values that
% change over time and do not have to restate all name-value pairs in
% .toNVpairs!
%
% After defining the animation struct, you can start the animation using
% the animateFromTo function:
%     animateFromTo(vidObj, animationStruct, duration, ...)
%
% See documentation for animateFromTo for more information.

%% Example 1: (simple)
% The following example animates a linear transformation of a zonotope 

% init set
rng(1)
S = zonotope.generateRandom('Dimension',3);
M = rand(3) * 4 -2;

% init animation
animationStruct = struct;
animationStruct(1).fromS = S;
animationStruct(1).dims = 1:3;
animationStruct(1).fromNVpairs = {'FaceColor',[1 0 0], 'FaceAlpha',0.2,'DisplayName','Z'};
animationStruct(1).transformFun = animateLinearTransform(M);
% animationStruct(1).toS = zonotope.generateRandom('Dimension',3); % alternative
animationStruct(1).toNVpairs = {'FaceColor',[0 0 1]};

animateFromTo(vidObj,animationStruct,2,'FromView',[-35,30],'ToView',[-15,20])

%% Example 2: (advanced)
% The following example animates an interval, an ellipsoid, and a zonotope,
% cycling through them and highlighting the currently selected set.
%
% The final video can be seen here: (TODO: update to youtube link)
%     https://drive.google.com/file/d/113QBER20Yzn9wVCWyS55b9diVpgDLzdK/view?usp=sharing

% interval
I = interval(zonotope([-3;0;0],eye(3)));

% ellipsoid
E = ellipsoid(eye(3),zeros(3,1));
E = zonotope(E,'outer:norm',28);

% zonotope
c = [ 3 ; 0 ; 0];
G = [ 0.123 0.018 0.080 0.556 -0.215 0.007 ; 0.568 -0.093 -0.054 -0.167 -0.028 -0.088 ; 0.190 0.122 -0.090 -0.008 0.237 -0.354 ];
Z = zonotope(c,G);

% name-value pairs for each set
NVpairsI = {'FaceColor',CORAcolor('CORA:color1'),'DisplayName','Interval'};
NVpairsE = {'FaceColor',CORAcolor('CORA:color2'),'DisplayName','Ellipsoid'};
NVpairsZ = {'FaceColor',CORAcolor('CORA:color3'),'DisplayName','Zonotope'};
NVpairsActive = {'FaceAlpha',0.25};
NVpairsInactive = {'FaceAlpha',0.05};
offSetActive = [0;0;0.5];

% define cycle
views = {[0,20],[35,30],[0,20],[-35,30],[0,20]};
sets = {I,E,Z};
NVpairs = {NVpairsI,NVpairsE,NVpairsZ};
idxActive =   [3,2,1,2];   % becomes active
idxInactive = [2,3,2,1];   % becomes inactive
idxPassive =  [1,1,3,3];   % stays passive

for i=1:numel(idxActive)
    % read respective sets
    setActive = sets{idxActive(i)};
    NVpairsSetActive = NVpairs{idxActive(i)};
    setInactive = sets{idxInactive(i)};
    NVpairsSetInactive = NVpairs{idxInactive(i)};
    setPassive = sets{idxPassive(i)};
    NVpairsPassive = NVpairs{idxPassive(i)};
    dims = 1:3;

    % build animation
    animationStruct = struct;
    % make active set active
    animationStruct(1).fromS = setActive;
    animationStruct(1).dims = dims;
    animationStruct(1).fromNVpairs = [NVpairsSetActive,NVpairsInactive];
    animationStruct(1).transformFun = animateLinearTransform(eye(3),offSetActive);
    animationStruct(1).toNVpairs = NVpairsActive;
    % make inactive set inactive
    animationStruct(2).fromS = setInactive + offSetActive;
    animationStruct(2).dims = dims;
    animationStruct(2).fromNVpairs = [NVpairsSetInactive,NVpairsActive];
    animationStruct(2).transformFun = animateLinearTransform(eye(3),-offSetActive);
    animationStruct(2).toNVpairs = NVpairsInactive;

    % plot passive set directly
    han = plot(setPassive,dims,NVpairsPassive{:},NVpairsInactive{:});

    % animate
    axis equal
    handles = animateFromTo(vidObj,animationStruct,2, ...
        "FromView",views{i},"ToView",views{i+1}, ...
        'FreezeDurationBefore',1,'FreezeDurationAfter',1, ...
        'FromEnlargeAxis',true,'ToEnlargeAxis',true);

    % delete all handles (removes all plotted sets for next animation)
    delete(han);
    for j=1:numel(handles)
        delete(handles{j});
    end
end

%% 3. FINISH CORA VIDEO ---------------------------------------------------

finishCORAvideo(vidObj,fig)

%% 4. (Optional) LOOP VIDEO AND ADD MUSIC ---------------------------------

% take a video, loop it, and add some background music
% a. install ffmpeg and make it available in the command line
% b. choose total duration of video below
% c. run this code snippet

% settings
filenameMusic = [CORAROOT '/global/functions/verbose/videos/July - John Patitucci.mp3'];
filenameVideoLooped = strrep(filenameVideo,'.mp4','_looped.mp4');
TotalDuration = 60;

% Read the video file
videoFReader = VideoReader(filenameVideo);
duration = videoFReader.Duration;

% determine number of loops
loops = ceil(TotalDuration / duration);

system( ...
    sprintf('ffmpeg -stream_loop %i -i "%s" -i "%s" -c:v copy -c:a aac -shortest %s', ...
    loops, filenameVideo, filenameMusic, filenameVideoLooped));

%% Stitch videos side by side

system(sprintf('ffmpeg -i ./CORAvideos/07/video.mp4 -i ./CORAvideos/07/example_angle.mp4 -filter_complex "hstack,format=yuv420p" -c:v libx264 -crf 18 video.mp4'))

%% 5. (Optional) SAVE TITLE, DESCRIPTION, AND THUMBNAIL -------------------

% readme ---

fileID = fopen([basepath '/readme.txt'],'w');

readmetext = "\n" + ...
    "Title:\n" + ...
    "[CORA] %s: %s\n" + ...
    "\n\n" + ...
    "Description:\n" + ...
    "%s\n" + ...
    "\n" + ...
    "For more information, please visit our website: https://cora.in.tum.de\n" + ...
    "\n\n" + ...
    "Keywords:\n" + ...
    "reachability analysis; formal verification; dynamic systems; differential equations; set-based computing; continuous sets\n" + ...
    "\n";

fprintf(fileID, readmetext, titleVideo, exampleName, compose(descVideo));
fclose(fileID);

% thumbnail ---

% Create a VideoReader object
v = VideoReader(filenameVideo);
FrameRate = v.FrameRate;

% Select frame for thumbnail (use Inf for last frame)
thumbnailFrame = read(v, Inf*FrameRate);

% Crop the last frame to desired dimensions (requires Image Processing Toolbox)
% Specify the rectangle for cropping as [x y width height]
croppedFrame = imcrop(thumbnailFrame, [2000 0 1527 1526]);

% Resize the cropped image if it's too large
maxWidth = 1280; % Maximum width for the thumbnail
maxHeight = 720; % Maximum height for the thumbnail

% Calculate the resize scale while keeping the aspect ratio
scale = min(maxWidth / size(croppedFrame, 2), maxHeight * 0.9 / size(croppedFrame, 1));
croppedFrame = imresize(croppedFrame, scale);

% Create a white canvas with the resolution of 1280x720 pixels
whiteCanvas = uint8(255 * ones(maxHeight, maxWidth, 3));

% Calculate the position to place the cropped image on the canvas
% This will center the cropped image on the canvas
startX = round((maxWidth - size(croppedFrame, 2)) / 2) +1;
startY = round((maxHeight - size(croppedFrame, 1)) / 2) +1;

% Place the cropped image onto the white canvas
whiteCanvas(startY:(startY + size(croppedFrame, 1) - 1), ...
startX:(startX + size(croppedFrame, 2) - 1), :) = croppedFrame;

% add TUM logo
logo = imread('TUM_logo.png','BackgroundColor',[1 1 1]);
logo = imresize(logo, 0.1);

startX = maxWidth - size(logo, 2)-25;
startY = 25;
whiteCanvas(startY:(startY + size(logo, 1) - 1), ...
startX:(startX + size(logo, 2) - 1), :) = logo;

% add CORA logo
logo = imread('coraLogo.png','BackgroundColor',[1 1 1]);
logo = imresize(logo, 0.175);
startX = 25;
startY = 25;
whiteCanvas(startY:(startY + size(logo, 1) - 1), ...
startX:(startX + size(logo, 2) - 1), :) = logo;

% Save the thumbnail to a file
imwrite(whiteCanvas, [basepath '/thumbnail.png']);

%% ------------------------------------------------------------------------

% ------------------------------ END OF CODE ------------------------------
