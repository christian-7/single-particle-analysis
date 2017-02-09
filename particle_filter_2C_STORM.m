% Input data from 2C particle segmentation

% Variable Cent

% For 2C Data

% 1 - locs              Channel 1
% 2 - number of locs in Channel 1
% 3 - int Intensity     Channel 1
% 4 - cropped WF ROI    Channel 1
% 5 - locs              Channel 2
% 6 - number of locs in Channel 2
% 7 - int Intensity     Channel 2
% 8 - cropped WF ROI    Channel 2


%% Filter by the number of locs per ROI

% Cent = all;

count = 0;
Cent_filt = {};

for i = 1:length(Cent);
    
    if isempty(Cent{i,1}) | isempty(Cent{i,5})
    
    else
        
NbrOfLocs(i,1) = length(Cent{i,1});

if length(Cent{i,1})>300 & length(Cent{i,1})<2e4;
   
count = count+1;    
    
Cent_filt{count,1} = Cent{i,1}; 
Cent_filt{count,2} = Cent{i,2};
Cent_filt{count,3} = Cent{i,3};
Cent_filt{count,4} = Cent{i,4};
Cent_filt{count,5} = Cent{i,5};
Cent_filt{count,6} = Cent{i,6};
Cent_filt{count,7} = Cent{i,7};
Cent_filt{count,8} = Cent{i,8};

else end

    end
end

figure
hist(NbrOfLocs,50);
xlabel('Number of localizations');
xlabel('Count');
title(['Median = ' num2str(median(NbrOfLocs))])

display([num2str(length(Cent_filt)/length(Cent)) '  left after filtering'])

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start FRC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute FRC resolution of a selected Particle

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cent_ID     = 6;  % Particle ID
Channel     = 1;    % Channel 1 - A647 / Channel 5 - A750
Channel2    = 5;

pixelsize   = 106;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coords_C1      = [];
coords_C1(:,1) = (Cent_filt{Cent_ID,Channel}(:,1) - min(Cent_filt{Cent_ID,Channel}(:,1)))/pixelsize;
coords_C1(:,2) = (Cent_filt{Cent_ID,Channel}(:,2) - min(Cent_filt{Cent_ID,Channel}(:,2)))/pixelsize;;
coords_C1(:,3) =  Cent_filt{Cent_ID,Channel}(:,4); % frame

coords_C2      = [];
coords_C2(:,1) = (Cent_filt{Cent_ID,Channel2}(:,1) - min(Cent_filt{Cent_ID,Channel2}(:,1)))/pixelsize;
coords_C2(:,2) = (Cent_filt{Cent_ID,Channel2}(:,2) - min(Cent_filt{Cent_ID,Channel2}(:,2)))/pixelsize;;
coords_C2(:,3) =  Cent_filt{Cent_ID,Channel2}(:,4); % frame


superzoom = 10;
szx_C1 = superzoom * round(max(coords_C1(:,1))-min(coords_C1(:,1)))*1.5;
szy_C1 = superzoom * round(max(coords_C1(:,2))-min(coords_C1(:,2)))*1.5;

szx_C2 = superzoom * round(max(coords_C2(:,1))-min(coords_C2(:,1)))*1.5;
szy_C2 = superzoom * round(max(coords_C2(:,2))-min(coords_C2(:,2)))*1.5;

im_C1 = binlocalizations(coords_C1, szx_C1, szy_C1, superzoom);
h=dipshow(im_C1);
dipmapping(h,[0 10],'colormap',hot)

im_C2 = binlocalizations(coords_C2, szx_C2, szy_C2, superzoom);
h=dipshow(im_C2);
dipmapping(h,[0 10],'colormap',hot)


%% compute resolution for both channels

fprintf('\n -- computing resolution Channel 1 --\n')
[res_value, ~, resH, resL] = postoresolution(coords_C1, szx_C1, superzoom); % in super-resolution pixels
fprintf('resolution value %2.1f +- %2.2f [px]\n', res_value, (resL-resH)/2);
fprintf('resolution value %2.1f +- %2.2f [nm]\n', res_value*pixelsize/superzoom, (resL-resH)/2*pixelsize/superzoom);

fprintf('\n -- computing resolution Channel 2 --\n')
[res_value, ~, resH, resL] = postoresolution(coords_C2, szx_C2, superzoom); % in super-resolution pixels
fprintf('resolution value %2.1f +- %2.2f [px]\n', res_value, (resL-resH)/2);
fprintf('resolution value %2.1f +- %2.2f [nm]\n', res_value*pixelsize/superzoom, (resL-resH)/2*pixelsize/superzoom);

% compute FRC curve
fprintf('\n -- computing FRC curve--\n')

[~,frc_curve_C1] = postoresolution(coords_C1, szx_C1, superzoom); 
[~,frc_curve_C2] = postoresolution(coords_C2, szx_C2, superzoom); 

figure('Position',[200 600 500 500]);hold all;
qmax = 0.5/(pixelsize/superzoom);
plot(linspace(0,qmax*sqrt(2), length(frc_curve_C1)), frc_curve_C1,'r-')
plot(linspace(0,qmax*sqrt(2), length(frc_curve_C2)), frc_curve_C2,'b-')
xlim([0,qmax])
hold on
plot([0 qmax],[1/7 1/7],'r-');
plot([0 qmax],[0 0],'k--');
xlabel('spatial frequency (nm^{-1})')
ylabel('FRC')
title('Fourier Ring Correlation curve')
legend('Channel 1', 'Channel 2');
box on

%% compute resolution (averaged over 20 runs)

Channel     = 2;    % Channel 1 - A647 / Channel 5 - A750

if Channel==1;

fprintf('\n -- computing resolution for Channel 1 averaged over 20 random block splits --\n')
[res_value, ~, resH, resL] = postoresolution(coords_C1, szx_C1, superzoom, 500,[], 20); % in super-resolution pixels
fprintf('res value %2.1f +- %2.2f [px]\n', res_value, (resL-resH)/2);
fprintf('res value %2.1f +- %2.2f [nm]\n', res_value*pixelsize/superzoom, (resL-resH)/2*pixelsize/superzoom);

else
    
fprintf('\n -- computing resolution for Channel 2 averaged over 20 random block splits --\n')
[res_value, ~, resH, resL] = postoresolution(coords_C2, szx_C2, superzoom, 500,[], 20); % in super-resolution pixels
fprintf('res value %2.1f +- %2.2f [px]\n', res_value, (resL-resH)/2);
fprintf('res value %2.1f +- %2.2f [nm]\n', res_value*pixelsize/superzoom, (resL-resH)/2*pixelsize/superzoom);
  
end

%%  compute resolution as a function of frame time (takes ~1-2min.)

Channel = 1;    % Channel 1 - A647 / Channel 5 - A750

if Channel==1;

fprintf('\n -- computing resolution for Channel 1 as a function of time for in 25 steps--\n')
tfrac = 25;
[~,~,~,~,resT] = postoresolution(coords_C1, szx_C1, superzoom, 500, tfrac); % in super-resolution pixels
figure;
plot(linspace(0,1,tfrac),resT*pixelsize/superzoom,'x-')
xlabel('time fraction')
ylabel('Resolution (nm)')
title('Resolution as a function of total number of frames')

else

fprintf('\n -- computing resolution for Channel 2 as a function of time for in 25 steps--\n')
tfrac = 25;
[~,~,~,~,resT] = postoresolution(coords_C2, szx_C2, superzoom, 500, tfrac); % in super-resolution pixels
figure;
plot(linspace(0,1,tfrac),resT*pixelsize/superzoom,'x-')
xlabel('time fraction')
ylabel('Resolution (nm)')
title('Resolution as a function of total number of frames')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End FRC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generate Variable for DBSCAN

Cent_filt_2C = {}; % New variable combining both datasets

for i = 1:length(Cent_filt)
    
Cent_filt{i,1}(:,10) = 1;
Cent_filt{i,5}(:,10) = 2;

Cent_filt_2C{i,1} = vertcat(Cent_filt{i,1},Cent_filt{i,5});

end

%% DBSCAN particle size filter

% % Output:   1. locs Ch1
%             2. locs Ch2
%             3. Rg Ch1
%             4. Rg Ch2
%             5. Ecc Ch1
%             6. Ecc Ch2

tic

DBSCAN_filtered = {};

for m = 1:length(Cent_filt_2C);
    
% [DBSCAN_res] = DBSCAN_batch(dataDBS,minLength);
% DBSCAN_res[locResults, Rg, Ecc]

[DBSCAN_filtered_temp] = DBSCAN_batch_2C(Cent_filt_2C{m,1},100);

DBSCAN_filtered = vertcat(DBSCAN_filtered,DBSCAN_filtered_temp);

clear DBSCAN_filtered_temp

clc
X = [' Finished DBSCAN ',num2str(m),' of ',num2str(length(Cent_filt_2C)),];
disp(X)

end

save('Test_2C_dataset.mat','DBSCAN_filtered');

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%%  Test DBSCAN particle size Filter

tic
clc, close all

m = 100; % Cent ID

fprintf('\n -- DBSCAN started --\n')
x = 1; y = 2;
DBSCAN_filtered = {};

% Find out if the data was processed in bstore or TS

dataDBS      = [];
dataDBS(:,1) = Cent_filt_2C{m,1}(:,1); % x in mum
dataDBS(:,2) = Cent_filt_2C{m,1}(:,2); % y in mum


if isempty(dataDBS)
   Cent_filt_2C{m,4} = [];
else

% Run DBSCAN on each particle 

k   = 10;                                                   % minimum number of neighbors within Eps
Eps = 15;                                                   % minimum distance between points, nm

[class,type]=DBSCAN(dataDBS,k,Eps);                         % uses parameters specified at input
class2=transpose(class);                                    % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                                      % (core: 1, border: 0, outlier: -1)

coreBorder = [];
coreBorder = find(type2 >= 0);

subset          = [];
subset          = Cent_filt_2C{m,1}(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);


subsetP = [];

subsetP(:,1)    = dataDBS(coreBorder,1);
subsetP(:,2)    = dataDBS(coreBorder,2);
subsetP(:,3)    = class2(coreBorder);

figure('Position',[700 600 1000 300])
subplot(1,3,1)
scatter(dataDBS(:,1),dataDBS(:,2),1);
title('Raw Data')
axis on
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

subplot(1,3,2)
scatter(subsetP(:,1),subsetP(:,2),1,mod(subsetP(:,3),10))
title('identified Clusters')
axis on
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

% Select only the largest cluster(s)

if isempty(subset);
else

ClusterLength = [];

for i = 1:max(subset(:,end));           % find the i-th cluster
    
    vx = find(subset(:,end)==i);
    
    if length(vx)>500;
        
    DBSCAN_filtered{length(DBSCAN_filtered)+1,1} = subset(vx,1:end);   
    
    subplot(1,3,3)
    scatter(DBSCAN_filtered{length(DBSCAN_filtered),1}(:,1),DBSCAN_filtered{length(DBSCAN_filtered),1}(:,2),1,'k');hold on;
    title('selected Clusters')
    axis on
    axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
    box on
    
    else end
   
end


end

end

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Compute FRC for all particles

pixelsize = 106; superzoom = 10; tic;

% Channel 1

for i = 1:length(DBSCAN_filtered);  

coords_C1      = [];
coords_C1(:,1) = (DBSCAN_filtered{i,1}(:,1) - min(DBSCAN_filtered{i,1}(:,1)))/pixelsize;
coords_C1(:,2) = (DBSCAN_filtered{i,1}(:,2) - min(DBSCAN_filtered{i,1}(:,2)))/pixelsize;
coords_C1(:,3) =  DBSCAN_filtered{i,1}(:,4) - min(DBSCAN_filtered{i,1}(:,4));

szx_C1 = superzoom * round(max(coords_C1(:,1))-min(coords_C1(:,1)))*1.5;
szy_C1 = superzoom * round(max(coords_C1(:,2))-min(coords_C1(:,2)))*1.5;

try
    
[res_value_C1, ~, resH_C1, resL_C1] = postoresolution(coords_C1, szx_C1, superzoom, 500,[], 5); % in super-resolution pixels

catch
    
 DBSCAN_filtered{i,7} = NaN;

continue
end

DBSCAN_filtered{i,7} = res_value_C1*pixelsize/superzoom; 

clc; 
X = [' Channel 1: Calculated ',num2str(i),' of ',num2str(length(DBSCAN_filtered)),];
clc;disp(X)

end

% Channel 2

for i = 1:length(DBSCAN_filtered);  

coords_C2      = [];
coords_C2(:,1) = (DBSCAN_filtered{i,2}(:,1) - min(DBSCAN_filtered{i,2}(:,1)))/pixelsize;
coords_C2(:,2) = (DBSCAN_filtered{i,2}(:,2) - min(DBSCAN_filtered{i,2}(:,2)))/pixelsize;
coords_C2(:,3) =  DBSCAN_filtered{i,2}(:,4) - min(DBSCAN_filtered{i,2}(:,4));

szx_C2 = superzoom * round(max(coords_C2(:,1))-min(coords_C2(:,1)))*1.5;
szy_C2 = superzoom * round(max(coords_C2(:,2))-min(coords_C2(:,2)))*1.5;

try
    
[res_value_C2, ~, resH_C2, resL_C2] = postoresolution(coords_C2, szx_C2, superzoom, 500,[], 5); % in super-resolution pixels

catch
    
    DBSCAN_filtered{i,8} = NaN;

continue
end

DBSCAN_filtered{i,8} = res_value_C2*pixelsize/superzoom; 

clc; 
X = [' Channel 2: Calculated ',num2str(i),' of ',num2str(length(DBSCAN_filtered)),];
clc;disp(X)

end

fprintf(' -- FRC computed in %f sec -- \n',toc)

save('2016-01-19_humanCent_aTub-Cep152_A647_all_extractedCent_filtered_DBSCAN.mat','DBSCAN_filtered');

%% Filter out NaN rows

Temp = [];
RowsToDelete = [];

for i = 1: length(DBSCAN_filtered)
    
    if or(isnan(DBSCAN_filtered{i,7}) == 1, isnan(DBSCAN_filtered{i,8})==1);
    
    RowsToDelete = vertcat(RowsToDelete,i);    
               
    elseif or(isempty(DBSCAN_filtered{i,7}) == 1, isempty(DBSCAN_filtered{i,8})==1);
    
    RowsToDelete = vertcat(RowsToDelete,i);     
        
    end
    
end

[Temp,ps] = removerows(DBSCAN_filtered,'ind',[RowsToDelete]);

DBSCAN_filtered = Temp;
clear Temp

toPrint = [' -- Cleaned DBSCAN filtered.  ' num2str(length(RowsToDelete)) ' Particles were filtered out --'];

disp(toPrint);

%% Compute Hollowness 

% DBSCAN_filtered=Cep152_all_filtered;

tic

for i = 1:length(DBSCAN_filtered);

[MeanH, StdH, MinH, MaxH] = calculate_Hollowness(DBSCAN_filtered{i,1});

DBSCAN_filtered{i,9} = MeanH;
DBSCAN_filtered{i,10} = StdH;
DBSCAN_filtered{i,11} = MinH;
DBSCAN_filtered{i,12} = MaxH;

end

for i = 1:length(DBSCAN_filtered);

[MeanH, StdH, MinH, MaxH] = calculate_Hollowness(DBSCAN_filtered{i,2});

DBSCAN_filtered{i,13} = MeanH;
DBSCAN_filtered{i,14} = StdH;
DBSCAN_filtered{i,15} = MinH;
DBSCAN_filtered{i,16} = MaxH;

end

clc
X = [' -- Calculated Hollowness descriptors in ', num2str(toc),' sec -- '];
disp(X)

%% Compute Circularity Ratio 

tic

for i = 1:length(DBSCAN_filtered);

[CircRatio] = calculate_Cicularity(DBSCAN_filtered{i,1});

DBSCAN_filtered{i,17} = CircRatio;

end

for i = 1:length(DBSCAN_filtered);

[CircRatio] = calculate_Cicularity(DBSCAN_filtered{i,2});

DBSCAN_filtered{i,18} = CircRatio;

end

clc

X = [' -- Calculated Circularity Ratio in ', num2str(toc),' sec -- '];
disp(X)

%% Put the data into a new nestes structure for clarity

for i = 1:length(DBSCAN_filtered);

Data_DBSCANed.Channel1.locs{i,1} = DBSCAN_filtered{i,1};
Data_DBSCANed.Channel2.locs{i,1} = DBSCAN_filtered{i,2};

Data_DBSCANed.Channel1.Rg{i,1} = DBSCAN_filtered{i,3};
Data_DBSCANed.Channel2.Rg{i,1} = DBSCAN_filtered{i,4};

Data_DBSCANed.Channel1.Ecc{i,1} = DBSCAN_filtered{i,5};
Data_DBSCANed.Channel2.Ecc{i,1} = DBSCAN_filtered{i,6};

Data_DBSCANed.Channel1.FRC{i,1} = DBSCAN_filtered{i,7};
Data_DBSCANed.Channel2.FRC{i,1} = DBSCAN_filtered{i,8};

Data_DBSCANed.Channel1.MeanHo{i,1} = DBSCAN_filtered{i,9};
Data_DBSCANed.Channel1.StdHo{i,1} = DBSCAN_filtered{i,10};
Data_DBSCANed.Channel1.MinHo{i,1} = DBSCAN_filtered{i,11};
Data_DBSCANed.Channel1.MaxHo{i,1} = DBSCAN_filtered{i,12};

Data_DBSCANed.Channel2.MeanHo{i,1} = DBSCAN_filtered{i,13};
Data_DBSCANed.Channel2.StdHo{i,1} = DBSCAN_filtered{i,14};
Data_DBSCANed.Channel2.MinHo{i,1} = DBSCAN_filtered{i,15};
Data_DBSCANed.Channel2.MaxHo{i,1} = DBSCAN_filtered{i,16};

Data_DBSCANed.Channel1.CircRatio{i,1} = DBSCAN_filtered{i,17};
Data_DBSCANed.Channel2.CircRatio{i,1} = DBSCAN_filtered{i,18};

end

save('2016-01-19_humanCent_aTub-Cep152_A647_all_extractedCent_filtered_DBSCAN.mat','Data_DBSCANed');

%% Plot Descriptors in subplots

close all

figure('Position',[100 100 1000 900],'Name','Channel1')

subplot(3,3,1)
hist(cell2mat(Data_DBSCANed.Channel1.Rg),30);
title('Radius of Gyration');
box on

subplot(3,3,2)
hist(cell2mat(Data_DBSCANed.Channel1.Ecc),30);
title('Eccentricity');
box on

subplot(3,3,3)
scatter(cell2mat(Data_DBSCANed.Channel1.Rg),cell2mat(Data_DBSCANed.Channel1.Ecc),3);
title('Rg vs Ecc');
xlabel('Radius of gyrration')
ylabel('Eccentricity')
box on

subplot(3,3,4)
hist(cell2mat(Data_DBSCANed.Channel1.FRC),30);
title(['FRC resolution, Median = ', num2str(median(cell2mat(Data_DBSCANed.Channel1.FRC)))]);
box on

subplot(3,3,5)
hist(cell2mat(Data_DBSCANed.Channel1.MeanHo),30);
title(['Mean Hollowness']);
box on

subplot(3,3,6)
hist(cell2mat(Data_DBSCANed.Channel1.StdHo),30);
title(['Std Hollowness']);
box on

subplot(3,3,7)
scatter(cell2mat(Data_DBSCANed.Channel1.MinHo),cell2mat(Data_DBSCANed.Channel1.MaxHo),3);
title('Min Hollown vs Max Hollown');
xlabel('Min Hollowness')
ylabel('Max Hollowness')
box on

subplot(3,3,8)
hist(cell2mat(Data_DBSCANed.Channel1.CircRatio),30);
title('Circularity');
box on

subplot(3,3,9)
scatter(cell2mat(Data_DBSCANed.Channel1.MeanHo),cell2mat(Data_DBSCANed.Channel1.StdHo),2);
title('Mean vs Std');
xlabel('Mean Hollowness')
ylabel('Std Hollowness')
box on

% Channel 2

figure('Position',[800 100 1000 900],'Name','Channel 2')

subplot(3,3,1)
hist(cell2mat(Data_DBSCANed.Channel2.Rg),30);
title('Radius of Gyration');
box on

subplot(3,3,2)
hist(cell2mat(Data_DBSCANed.Channel2.Ecc),30);
title('Eccentricity');
box on

subplot(3,3,3)
scatter(cell2mat(Data_DBSCANed.Channel2.Rg),cell2mat(Data_DBSCANed.Channel2.Ecc),3);
title('Rg vs Ecc');
xlabel('Radius of gyrration')
ylabel('Eccentricity')
box on

subplot(3,3,4)
hist(cell2mat(Data_DBSCANed.Channel2.FRC),30);
title(['FRC resolution, Median = ', num2str(median(cell2mat(Data_DBSCANed.Channel2.FRC)))]);
box on

subplot(3,3,5)
hist(cell2mat(Data_DBSCANed.Channel2.MeanHo),30);
title(['Mean Hollowness']);
box on

subplot(3,3,6)
hist(cell2mat(Data_DBSCANed.Channel2.StdHo),30);
title(['Std Hollowness']);
box on

subplot(3,3,7)
scatter(cell2mat(Data_DBSCANed.Channel2.MinHo),cell2mat(Data_DBSCANed.Channel2.MaxHo),3);
title('Min Hollown vs Max Hollown');
xlabel('Min Hollowness')
ylabel('Max Hollowness')
box on

subplot(3,3,8)
hist(cell2mat(Data_DBSCANed.Channel2.CircRatio),30);
title('Circularity');
box on

subplot(3,3,9)
scatter(cell2mat(Data_DBSCANed.Channel2.MeanHo),cell2mat(Data_DBSCANed.Channel2.StdHo),2);
title('Mean vs Std');
xlabel('Mean Hollowness')
ylabel('Std Hollowness')
box on


%% Filter by length, Rg, Ecc and FRC

count = 0;
Cent_selected = {};

MinLength = 200;
MinRg = 50; MaxRg = 200;
MinEcc = 1; MaxEcc = 3;
MinFRC = 100;
% First filter by length

for i = 1:length(DBSCAN_filtered);
    
if length(DBSCAN_filtered{i,1})>MinLength   == 1        & ... % length Ch1
   length(DBSCAN_filtered{i,2})>MinLength   == 1        & ... % length Ch2
   DBSCAN_filtered{i,3}>MinRg               == 1        & ... % Rg Ch1
   DBSCAN_filtered{i,3}<MaxRg               == 1        & ... % Rg Ch1
   DBSCAN_filtered{i,4}>MinRg               == 1        & ... % Rg Ch2
   DBSCAN_filtered{i,4}<MaxRg               == 1        & ... % Rg Ch2
   DBSCAN_filtered{i,5}>MinEcc              == 1        & ... % Ecc Ch1
   DBSCAN_filtered{i,5}<MaxEcc              == 1        & ... % Ecc Ch1
   DBSCAN_filtered{i,6}>MinEcc              == 1        & ... % Ecc Ch2
   DBSCAN_filtered{i,6}<MaxEcc              == 1        & ... % Ecc Ch2
   DBSCAN_filtered{i,7}<MinFRC              == 1        & ... % Ecc Ch2
   DBSCAN_filtered{i,8}<MinFRC              == 1;
   
count = count+1;    
    
Cent_selected{count,1} = DBSCAN_filtered{i,1}; 
Cent_selected{count,2} = DBSCAN_filtered{i,2};

else end

end

X = [num2str(length(Cent_selected)) '/' num2str(length(DBSCAN_filtered)) '  particles are left'];
disp(X)

%% Render the filtered images 

pxlsize = 10;

save_path = 'Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\2C STORM analysis\rendered_images';
cd(save_path);

% Determine the box size form the largest particle

im_size = []; x = 1; y = 2;

for i=1:length(Cent_selected);
           
im_size(i,1) = round((max(Cent_selected{i,1}(:,y))-min(Cent_selected{i,1}(:,y)))/pxlsize);
im_size(i,2) = round((max(Cent_selected{i,1}(:,x))-min(Cent_selected{i,1}(:,x)))/pxlsize);

im_size_C2(i,1) = round((max(Cent_selected{i,2}(:,y))-min(Cent_selected{i,2}(:,y)))/pxlsize);
im_size_C2(i,2) = round((max(Cent_selected{i,2}(:,x))-min(Cent_selected{i,2}(:,x)))/pxlsize);
              
end

% Images Channel 1

count=1;

for i = 1:length(Cent_selected);
    
        heigth = round((max(Cent_selected{i,1}(:,y)) - min(Cent_selected{i,1}(:,y)))/pxlsize);
        width  = round((max(Cent_selected{i,1}(:,x)) - min(Cent_selected{i,1}(:,x)))/pxlsize);
        
        rendered = hist3([Cent_selected{i,1}(:,y),Cent_selected{i,1}(:,x)],[heigth width]);
        empty  = zeros(round(max(max(im_size_C2))*1.5),round(max(max(im_size_C2))*1.5));

        size_DC    = size(empty);
        center_X   = round(size_DC(2)/2);
        center_Y   = round(size_DC(1)/2);

        empty(round(center_Y-heigth/2):(round(center_Y-heigth/2))+heigth-1,round(center_X-width/2):(round(center_X-width/2))+width-1) = rendered;

        name32rendered  = ['image_particles_Channel_1_PxlSize_' num2str(pxlsize) '_' num2str(i) '.tiff'];


I32 = [];
I32 = uint32(empty);

t = Tiff(name32rendered,'w');
tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()

count=count+1;

end

% Images Channel 2

count=1;

for i = 1:length(Cent_selected);
    
        heigth = round((max(Cent_selected{i,2}(:,y)) - min(Cent_selected{i,2}(:,y)))/pxlsize);
        width  = round((max(Cent_selected{i,2}(:,x)) - min(Cent_selected{i,2}(:,x)))/pxlsize);
        
        rendered = hist3([Cent_selected{i,2}(:,y),Cent_selected{i,2}(:,x)],[heigth width]);
        empty  = zeros(round(max(max(im_size_C2))*1.5),round(max(max(im_size_C2))*1.5));

        size_DC    = size(empty);
        center_X   = round(size_DC(2)/2);
        center_Y   = round(size_DC(1)/2);

        empty(round(center_Y-heigth/2):(round(center_Y-heigth/2))+heigth-1,round(center_X-width/2):(round(center_X-width/2))+width-1) = rendered;

        name32rendered  = ['image_particles_Channel_2_PxlSize_' num2str(pxlsize) '_' num2str(i) '.tiff'];


I32 = [];
I32 = uint32(empty);

t = Tiff(name32rendered,'w');
tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()

count=count+1;

end

clc;
disp([' -- Images Saved --']);

%% Use ImageJ to create, blur and save the Montage

% addpath('D:\fiji-win64\Fiji.app\scripts');
% Miji;

cd(save_path);
string1 = ['Image Sequence...']; string2 = ['open=' num2str(save_path) ' file=Channel1 sort'];

MIJ.run(string1, string2);
MIJ.run('Gaussian Blur...','sigma=1 stack');

string1 = ['Make Montage...']; string2 = ['columns=' num2str(round(sqrt(length(Cent_selected)))) ' rows=' num2str(round(sqrt(length(Cent_selected))))  ' scale=1'];

MIJ.run(string1, string2);
MIJ.run('Enhance Contrast', 'saturated=0.35');
MIJ.run('Save', 'Tiff..., path=[Z:\\Christian-Sieben\\data_HTP\\2016-04-01_humanCentriole_aTubNB_Sas6\\2C STORM analysis\\Montage_Channel1.tif]');
MIJ.run('Close All')

%% First Classification - TopOrSide / Unclassified

% This classifier uses the function waitforbuttonpress

% Start this after the density filtering

Cent_selected = DBSCAN_filtered(1:300,1:end);

Cent = Cent_selected;

Top                 = {};
unclassified        = {};
response = [];
          
% Top = {};
% Side= {};

clc, close all

for i = 1:length(Cent);

            
        figure('Position',[500 300 300 300]);
        
        scatter(Cent{i,1}(:,1)-min(Cent{i,1}(:,1)),Cent{i,1}(:,2)-min(Cent{i,1}(:,2)),1,'filled','black');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);
        axis([0 max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)) 0 max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2))])
        box on
        axis on  
        
        fprintf('Mouse --> Top or Side');
        fprintf('\n');
        fprintf('\n');
        fprintf('Key --> unclassified');
                
        w = waitforbuttonpress;                                     % 0 if it detects a mouse button click / 1 if it detects a key press
        
        if w == 0;  
            
        Top = vertcat(Top,Cent{i,1});              % Put on the shown cluster into the new variable  
        response(i,1) = 1; 
        disp('Top')
        
        else
            
        unclassified = vertcat(unclassified,Cent{i,1});        % Put on the shown cluster into the new variable 
        response(i,1) = 2;
        disp('Side')
        end


    close all  , clc     

       
end
%% First Classification - Top / Side / Unclassified

% This classifier uses the function input
% Start this after the density filtering

% Cent_selected = Cent_filt(1:300,1:end);

Cent_selected = DBSCAN_filtered(1:300,1:end);

Cent = Cent_selected;

Top             = {};
Side            = {};
Unclassified    = {};
response        = []; % 1>Top, 2>Side, 3>Unclassified

clc, close all

for i = 1:length(Cent);

            
        figure('Position',[500 300 300 300]);
        
        scatter(Cent{i,1}(:,1)-min(Cent{i,1}(:,1)),Cent{i,1}(:,2)-min(Cent{i,1}(:,2)),1,'filled','black');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);
        axis([0 max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)) 0 max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2))])
        box on
        axis on  

        w = input('Waiting for Input: t > Top, s > Side, u > unclassified','s')
                
        if strcmpi(w,'t');  
            
        Top = vertcat(Top,Cent{i,1});              % Put on the shown cluster into the new variable  
        response(i,1) = 1; 
        disp('Top View')
        
        elseif strcmpi(w,'s')
            
        Side = vertcat(Side,Cent{i,1});        % Put on the shown cluster into the new variable 
        response(i,1) = 2;
        disp('Side View')
        
        else 
            
        Unclassified = vertcat(Unclassified,Cent{i,1});        % Put on the shown cluster into the new variable 
        response(i,1) = 3;
        disp('Unclassified')
        
        end


    close all  , clc     

       
end




