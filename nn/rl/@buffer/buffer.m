classdef buffer
% buffer - replay buffer of the rl agent.
%   the array contains (s_i,a_i,r_i,s_i+1) and visualData for
%   plotting of trajectories during training process
%
% Syntax:
%   obj = buffer(bufferSize)
%
% Inputs:
%   bufferSize - maximal size of buffer
% 
% Outputs:
%   obj - instantiated buffer object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DDPGagent

% Authors:       Manuel Wendl
% Written:       24-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    properties
        array
        visualisationData
        bufferSize
        currentIndx
    end
    
    methods
        % constructor 
        function obj = buffer(bufferSize)
            obj.array = {};
            obj.visualisationData.episodeNum = [];
            obj.visualisationData.reachSet = {};
            obj.bufferSize = bufferSize;
            obj.currentIndx = 1;
        end

        function obj = resetBuffer(obj)
            obj.array = {[],[],[],[],[]};
            obj.visualisationData.episodeNum = [];
            obj.visualisationData.reachSet = {};
            obj.currentIndx = 1;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
