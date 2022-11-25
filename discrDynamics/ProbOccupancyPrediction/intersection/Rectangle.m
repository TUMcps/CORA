classdef Rectangle < handle
    %RECTANGLE Summary of this class goes here
    %   Detailed explanation goes here
    %width along Oy
    %height along Ox
    %4----1
    %|    |
    %3----2
    %update: 19-August-2016, Matthias Althoff
    
    properties (SetAccess = private)
        width
        height
        orientation
        position
        
    end
    properties (SetAccess = private, GetAccess = private)
        vertices
        radius
    end
    properties
        color = [0 0 1] % only for visualization
    end
    
    methods
        function obj = Rectangle(width, height, orientation, position)
            obj.width = width;
            obj.height = height;
            obj.orientation = orientation;
            obj.position = position;
            obj.radius = sqrt(width^2/4+height^2/4);
            
            w = obj.width/2;
            h = obj.height/2;
            vert = [w w -w -w; h -h -h h];
            vert = rot2d(obj.orientation, vert);
            vert = [vert(1,:) + obj.position(1); vert(2,:) + obj.position(2)];
            obj.vertices = vert;
        end
        
        function v = getver(obj)
            v = obj.vertices;
        end
        
        function r = getradius(obj)
            r = obj.radius;
        end
        
        % uses separating axis theorem
        function i = intersect(obj, rect)
            i = true;
            
            % check distance between center of both rectangles
%             if norm(rect.position-obj.position) >= rect.radius+obj.radius
%                 i = false;
%                 return;
%             end
            
            % initialize four candidate axis
            ax_l = [cos(obj.orientation) cos(obj.orientation+pi/2) ...
                cos(rect.orientation) cos(rect.orientation+pi/2); ...
                sin(obj.orientation) sin(obj.orientation+pi/2) ...
                sin(rect.orientation) sin(rect.orientation+pi/2)];
            
            % this is a bit faster
            %             ax_l = [obj.vertices(:,1)-obj.vertices(:,2), ...
            %                 obj.vertices(:,2)-obj.vertices(:,3), ...
            %                 rect.vertices(:,1)-rect.vertices(:,2), ...
            %                 rect.vertices(:,2)-rect.vertices(:,3)];
            for k=1:4
                ax = ax_l(:,k);
                p_sa = (ax'*obj.vertices); % project points of rect "a" on axis
                p_sb = (ax'*rect.vertices); % project points of rect "b" on axis
                s_a = sort(p_sa);
                s_b = sort(p_sb);
                
                % check if points are separated
                if s_a(1)>s_b(end) || s_a(end)<s_b(1) 
                    i = false;
                    return;
                end
            end
        end
               
        function draw(obj, type)
            v = obj.vertices;
            v = [v v(:,1)];
            plot(v(1,:), v(2,:), type);
        end
        
    end
    
end
