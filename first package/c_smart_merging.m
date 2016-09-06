function [Cent] = c_smart_merging(Locpath,locName,Segpath,segName);

% Input: segmented centrioles
% 
% Uses the segmeneted centrioles and performs smart merging on each segment
% 
% Output: segmented and merged centrioles --> It overwrites the .mat file,
% while adding a column 

%% Load Data and Header

% Load the header

cd(Locpath);

file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('x [nm]',h);
yCol = strmatch('y [nm]',h);
frameCol = strmatch('frame',h);
photonsCol = strmatch('intensity [photon]',h);
sigmaCol = strmatch('sigma [nm]',h);
uncertaintyCol = strmatch('uncertainty [nm]',h);

% Load the segemented centrioles

cd(Segpath);
load(segName);

% Load the fitting results

load('Z:\Christian-Sieben\software\smart_merging\fit_gap_scatter.mat')

fprintf('-- 1/6 Header and Data Loaded --\n')

for k = 1:length(Cent);
    
    sizeCent = size(Cent{k,1});
    
    if sizeCent(:,1) < 5;       % Tracker doesnt work for short segments
        Cent{k,2} = Cent{k,1};
    else
     
%% Create input fo the tracker

pos_list = [];

pos_list(:,1)=Cent{k,1}(:,xCol);                   % in pxl
pos_list(:,2)=Cent{k,1}(:,yCol);                   % in pxl
pos_list(:,3)=Cent{k,1}(:,photonsCol);             % photons
pos_list(:,4)=Cent{k,1}(:,uncertaintyCol);         % uncertainty
pos_list(:,5)=Cent{k,1}(:,frameCol);               % dt in seconds

fprintf('-- 2/6 Pos_list generated --\n')

%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

max_disp    = 15;           % in unit of data
min_pos     = 1;            % good - eliminate if fewer than good valid positions
gap         = 100;          % mem - number of time steps that a particle can be 'lost' and then recovered again
quiet       = 1;            % quiet - 1 = no text

param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res=trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID

fprintf('-- 3/6 Tracking Done --\n')

%% Calculate the distance and the gap time from origin of track

res(:,7)=zeros(1); % distance r from the center of the track
res(:,8)=zeros(1); % gap time between 1st and n-th localization

for i=1:max(res(:,6));          % for all tracks
    
    vx=find(res(:,6)==i);       % find the i-th track 

    track=res(vx,1:6);          % single track
    
    track_center(:,1) = sum(track(:,1))/length(track(:,1)); % center of the i-th track, x
    track_center(:,2) = sum(track(:,2))/length(track(:,2)); % center of the i-th track, y
    
    track(:,7)=zeros(1);
    track(:,8)=zeros(1);
    
    if length(track(:,1))>1;
    
    for j=1:length(track(:,1));
        
        track(j,7)= sqrt(((track_center(:,1)-track(j,1))^2)+((track_center(:,2)-track(j,2))^2)); % calculate distance from track_centers
        track(j,8)= track(j,5)-track(1,5); % time gap in frames
        
    end
    
    else end
    
    res(vx,7)=track(:,7); % distance
    res(vx,8)=track(:,8); % gap
  
end

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID
% res7 = distance
% res8 = gap time

fprintf('-- 4/6 Gap distance and time calculated --\n')

%% Assign Probability to each localization per cluster

res(:,9)=zeros(1);      % total probability = probablity from distance * probablity from time

sizeRes = size(res);

for i=1:sizeRes(:,1);  
   
    res(i,9)=fitresult(res(i,7))*gapfit(res(i,8)); % 

end

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID
% res7 = distance
% res8 = gap time
% res9 = probability

fprintf('-- 5/6 Probability assigned --\n')

%% Weighted Merging according to the probability

% Assign all empty variables

    x = [];
    y = [];
    Averagex = [];
    Averagey = [];
    
    groupedx = [];
    groupedy = [];
    groupedframe = [];
    groupedID = [];
    groupedPhotons = [];
    merged = [];
    groupedUncertainty = [];
    
% Perform merging according to probability
    

for i=1:max(res(:,6));          % for all tracks
    
    vx=find(res(:,6)==i);       % find the i-th track 

    track=res(vx,1:end);        % isolate single track --> track
    
    ProbFrac = (track(:,9)/(sum(track(:,9))))*100;  % calculate probability fraction (frac of 100)
    
    for j=1:length(vx);
       
    x(j) = track(j,1)*ProbFrac(j);  % x coordinate
    y(j) = track(j,2)*ProbFrac(j);  % y coordinate
    Averagex = sum(x)/100;            % calculate the average
    Averagey = sum(y)/100;            % calculate the average
    
    end
    
            groupedx=vertcat(groupedx,Averagex);
            groupedy=vertcat(groupedy,Averagey);
            groupedframe=vertcat(groupedframe, round(mean(track(:,5))));
            groupedID=vertcat(groupedID, i);
            groupedPhotons=vertcat(groupedPhotons, sum(track(:,3)));
            groupedUncertainty=vertcat(groupedUncertainty, mean(track(:,4)));
            
    x = [];
    y = [];
    Averagex = [];
    Averagey = [];
    
    
end

% x,y,photons,uncertainty, frame, ID

Cent{k,2}(:,1) = groupedx;
Cent{k,2}(:,2) = groupedy;
Cent{k,2}(:,3) = groupedPhotons;
Cent{k,2}(:,4) = groupedUncertainty;
Cent{k,2}(:,5) = groupedframe;
Cent{k,2}(:,6) = groupedID;

fprintf('-- 6/6 Molecules Merged  --\n');

fprintf('\n -- Finished Merging %f of %f  --\n', k,length(Cent));

    end
end

cd(Segpath);
save(segName,'Cent');

fprintf('\n -- File Saved --\n');
