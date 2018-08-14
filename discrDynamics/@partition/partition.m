classdef partition
%% class definition of a partition
% a partition divides the state-space, which is an axis-aligned bounding 
% box, into segments which are axis-aligned bounding boxes.
% The partition can be defined by one argument: a cell array of the values
% of the dividers in each dimension, or by two arguments: the state-space
% and the number of segments in each direction. For the latter definition,
% the dividers are evenly spaced along each dimension.
% 
% Usage:
% C = partition(dividers) -- where dividers is a cell array of arrays which
% specify the boundaries of the segments in each dimension
% C = partition(intervals, nrOfSegments) -- intervals is an n x 2 matrix of
% the state space, or an interval vector. nrOfSegments is an n x 1 vector
% of the number of partition segments along each dimension.
%
% See unit tests for examples.
%
% Last modified: AP 2.8.17




    properties
        intervals           % state space
        nrOfSegments        % formerly maxIndex
        % actualSegmentNr     % deprecated -- we handle possibly more than
        % 1 current segment now and variable-size fields are not permitted
        % in real-time Simulink.
        % mode                % deprecated
        dividers            % dividers
    end
    
    methods
        
        function Obj = partition(varargin)
            % empty object
            if nargin == 0
                disp('Partition needs more input values');

            % Else if the parameter is an identical object, copy object    
            elseif isa(varargin{1}, 'partition')
                Obj = varargin{1};
                
            % If two arguments are passed 
            elseif nargin == 2
                %=======================================================
                Obj.nrOfSegments=varargin{2};
                
                if isa(varargin{1},'interval')
                    if size(varargin{1},2) == 1     % if it is a column vector
                        Obj.intervals = [infimum(varargin{1}),supremum(varargin{1})];
                    else                            % make it a column vector!
                        Obj.intervals = [infimum(varargin{1})',supremum(varargin{1})'];
                    end
                else
                    Obj.intervals=varargin{1};
                end
                
                if size(Obj.nrOfSegments,1)~=size(Obj.intervals,1)
                    disp('invalid partition -- dimension mismatch or number of segments given as a row vector (should be column)!')
                    Obj.intervals = [];
                    Obj.nrOfSegments = [];
                end
                
                for i = 1:length(Obj.nrOfSegments)
                    Obj.dividers{i} = linspace(Obj.intervals(i,1),Obj.intervals(i,2),Obj.nrOfSegments(i)+1);
                end

            % If one argument is passed 
            elseif nargin == 1
                %=======================================================
                Obj.dividers=varargin{1};
                sup=cellfun(@(v) v(end), varargin{1}(:));
                inf=cellfun(@(v) v(1), varargin{1}(:));
                Obj.intervals = [inf,sup];
                for i=1:length(Obj.dividers)
                    Obj.nrOfSegments(i,1)=length(Obj.dividers{i})-1;
                end
                %=======================================================


            % Otherwise use a specific constructor    
            else
                disp('Partition needs more/fewer input values');
                Obj=[];
            end
            
            if any(Obj.nrOfSegments)<1
                disp('empty partition!')
            end
        end
        
        segmentMatrix = cellSegments(obj,indexVector)
        [segments,error]=intersectingCells(Obj,contSet,varargin)
        [cellIndices]=cellIndices(Obj,subscriptVectors)
        [c]=cellCenter(Obj,cellNr)
        intervals = cellIntervals(obj,varargin)
        polytopes = cellPolytopes(obj,varargin)
        zonotopes = cellZonotopes(obj,varargin)
        [intersectingCells, intersectingProportion]=exactIntersectingCells(Obj,contSet)
        n = nrOfCells(obj)
        [c]=centerSegment(Obj)
        [normalizedVector]=normalize(Obj,vector)
        plot(Obj,varargin)
        display(obj)
    end

end