function plot(varargin)
% Purpose:  plot road object
% Pre:      road object
% Post:     ---
% Built:    04.12.06,MA
% Modified: 21.11.07,MA
% Modified: 13.08.18,MA


%no color specified
if nargin==2
    obj=varargin{1}; %road object
    p=varargin{2}; %segment probabilities
    devProb=[1,2,3,3.5,3,2,1]/sum([1,2,3,3.5,3,2,1]); %deviation probability  
    color='b'; %b=blue
    %get maximum probability for normalization
    pMax=max(p)*1.02*max(devProb);
    normVal=1/pMax;     
    segment=[];
    
%deviation probability specified
elseif nargin==3
    obj=varargin{1};
    p=varargin{2}; %segment probabilities
    devProb=varargin{3};
    color='b'; %b=blue
    %get maximum probability for normalization
    pMax=max(p)*1.02*max(devProb);
    normVal=1/pMax;     
    segment=[];    
    
%color specified
elseif nargin==4
    obj=varargin{1};
    p=varargin{2}; %segment probabilities
    devProb=varargin{3};
    color=varargin{4};
    %get maximum probability for normalization
    pMax=max(p)*1.02*max(devProb);
    normVal=1/pMax;    
    segment=[];
    
%color, normVal specified
elseif nargin==5
    obj=varargin{1};
    p=varargin{2}; %segment probabilities
    devProb=varargin{3};
    color=varargin{4};
    normVal=varargin{5};
    segment=[];       
    
%color, normVal and segmentn specified
elseif nargin==6
    obj=varargin{1};
    p=varargin{2}; %segment probabilities
    devProb=varargin{3};
    color=varargin{4};
    normVal=varargin{5};
    segment=varargin{6};    
end


%get number of deviation segments
nrOfDev=obj.nrOfDevSegments;

%plot(x,y,'r+');

%find segments that should be plotted
if isempty(segment)
    ind=find(p);
else
    ind=segment;
end

for iInd=1:length(ind)
    
    %get segment number
    iSeg=ind(iInd);

    if p(iSeg)>0.0001
        %iDev: deviation segment
        for iDev=1:nrOfDev

            %obtain segment polytope
            [P]=segPolytope(obj,iSeg,iDev);
                        
            if strcmp(color,'trans')
                options.color='k';               
            elseif strcmp(color,'transB')
                options.color='b';
            elseif strcmp(color,'transG')
                options.color='g';
            elseif strcmp(color,'transR')
                options.color='r'; 
            else
                options.color=color;  
            end
            
            % no borders
            options.linestyle='none';
            
            % obtain transparency
            options.shade=p(iSeg)*normVal*devProb(iDev);
            
            % make small probabilities visible
            if options.shade<1e-6
                options.shade=1e-6;
            end
            
            % plot polytope
            try %MPT2
                plot(P,options);
            catch %MPT3
                try
                plot(P,'Color',options.color,'LineStyle',options.linestyle,'Alpha',options.shade);
                catch
                    disp('display error')
                end
            end

            hold on
        end
    end
end
