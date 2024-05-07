%% HOW TO RECORD VIDEOS IN CORA
%
% 1. Gather reachable sets and simulations
%   (e.g. set breakpoint after everything is computed / at the end)
% 2. Choose title and description (adapt example description)
% 3. Record video
% 4. Loop video and add music
%
% The code snippets below should give you a good starting point to create
%   your own video
%
% See also: recordCORAvideo

% Authors:       Tobias Ladner
% Written:       12-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% 1. GATHER REACHABLE SETS AND SIMULATIONS

% (e.g. set breakpoint after everything is computed / at the end)

% settings
basepath = './CORAvideos/0x';
exampleName = 'Van-der-Pol Oscillator';
mkdir(basepath)
filenameVideo = [basepath '/video.mp4'];

%% 2. TITLE AND DESCRIPTION -----------------------------------------------

% set title and description

%% Reachability Analysis for Continuous Systems

title = 'Reachability Analysis for Continuous Systems';

desc = "CORA computes reachable sets for linear systems, nonlinear systems as well as for systems with constraints. " + ...
    "Continuous as well as discrete time models are supported. Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitly considered. " + ...
    "In addition, CORA also provides capabilities for the simulation of dynamical models.\n" + ...
    "\n" + ...
    "The example on the right formally verifies the stability of the Van-der-Pol oscillator system. " + ...
    "For details of the system, please visit Sec. 3.1 of the ARCH-COMP19 category report: Continuous and hybrid systems with nonlinear dynamics.";


%% Reachability Analysis for Hybrid Systems

title = 'Reachability Analysis for Hybrid Systems';

desc = "CORA is capable to calculate the reachable sets for hybrid systems. " + ...
    "All implemented dynamic system classes can be used to describe the different continuous flows for the discrete system states. " + ...
    "Further, multiple different methods for the calculation of the intersections with guard sets are implemented in CORA.\n" + ...
    "\n" + ...
    "The example on the right shows a ball bouncing off a table over time. " + ...
    "This benchmark is particularly difficult to verify due to the instant transitions every time the ball hits the table.";


%% Formal Verification of Neural Networks

title = 'Formal Verification of Neural Networks';

desc = "CORA enables the formal verification of neural networks, both in open-loop as well as in closed loop scenarios. " + ...
    "Open-loop verification refers to the task where properties of the output set of a neural network are verified, e.g. correctly classified images given noisy input. " + ...
    "In closed-loop scenarios, the neural network is used as controller of a dynamic system and is neatly integrated in the reachability algorithms above, e.g. controlling a car while keeping a safe distance.\n" + ...
    "\n" + ...
    "The example on the right shows a car controlled by a neural network following a leading car. " + ...
    "We can formally verify that the NN-controlled car maintains a safe distance to the leading car.";


%% 3. RECORD VIDEO --------------------------------------------------------

% record video

%% plot

dims = 1:2;
recordCORAvideo(filenameVideo, ...
    'ReachSet',R,'SimResult',simRes,'Dimensions',dims, ...
    'Title', title,'Description', desc, 'Specification', spec)

%% plotOverTime

dims = 1;
recordCORAvideo(filenameVideo, ...
    'ReachSet',R,'SimResult',simRes,'Dimensions',dims, ...
    'Title', title,'Description', desc, ...
    'XLim', [0,params.tFinal], 'YLim', [0,1.2])


%% ACC

dims = 1;
recordCORAvideo(filenameVideo, ...
    'ReachSets',{ ...
        project(R_distances,1),'Distance',CORAcolor('CORA:safe'); ...
        project(R_distances,2),'Safe distance',CORAcolor('CORA:unsafe') ...
    },'SimResult',simResDistances,'Dimensions', dims, ...
    'Title', title,'Description', desc, ...
    'YLim', [0,110]);

%% 4. LOOP VIDEO AND ADD MUSIC --------------------------------------------

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

%%

% stitch videos side by side

system(sprintf('ffmpeg -i ./CORAvideos/07/video.mp4 -i ./CORAvideos/07/example_angle.mp4 -filter_complex "hstack,format=yuv420p" -c:v libx264 -crf 18 video.mp4'))

%% %. (Optional) SAVE TITLE, DESCRIPTION, and THUMBNAIL

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

fprintf(fileID, readmetext, title, exampleName, compose(desc));
fclose(fileID);

% thumbnail ---

% Create a VideoReader object
v = VideoReader(filenameVideo);

% Read the last frame
lastFrame = read(v, inf);

% Crop the last frame to desired dimensions
% Specify the rectangle for cropping as [x y width height]
croppedFrame = imcrop(lastFrame, [2000 0 1527 1526]);

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

% add TUM cps logo
logo = imread('TUM_logo.png','BackgroundColor',[1 1 1]);
logo = imresize(logo, 0.1);

startX = maxWidth - size(logo, 2)-25;
startY = 25;
whiteCanvas(startY:(startY + size(logo, 1) - 1), ...
startX:(startX + size(logo, 2) - 1), :) = logo;

% add cora logo
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
