%% Filter by the number of locs per ROI

Cent = all_part;

Cent_filt = {};

count = 0;

for i = 1:length(Cent);

if length(Cent{i,1})>300;
    
count = count+1;    
    
Cent_filt{count,1} = Cent{i,1};
Cent_filt{count,2} = Cent{i,2};
Cent_filt{count,3} = Cent{i,3};
Cent_filt{count,4} = Cent{i,4};
Cent_filt{count,5} = Cent{i,5};
Cent_filt{count,6} = Cent{i,6};

else end

end


%% NN filter

x=1; y=2;

for j = 1;
    
    NN = rangesearch(Cent_filt{j,1}(:,x:y),Cent_filt{j,1}(:,x:y),50);
    
end

NofN = [];

for i = 1:length(NN);
    
    NofN(i,1) = length(NN{i,1});

end

minNN = 50;

figure
scatter(Cent_filt{j,1}(:,x)/1000,Cent_filt{j,1}(:,y)/1000,1,'black');hold on;
scatter(Cent_filt{j,1}(NofN>minNN,x)/1000,Cent_filt{j,1}(NofN>minNN,y)/1000,1,'o','red');

%%  DBSCAN particle size Filter

tic
fprintf('\n -- DBSCAN started --\n')

% Delete the 6th column of the Structure

for i=1:length(Cent_filt);
    
Cent_filt{i,6} = [];

end

% Find out if the data was processed in bstore or TS

for m = 1:length(Cent_filt);

if max(Cent_filt{m,1}(:,4)) == 0 %max(Cent_filt{m,1}(:,1))/max(Cent_filt{m,1}(:,2))>1    
    
    x = 2; y = 3;

else
    
    x = 1; y = 2;

end

% Select the data for DBSCAN

dataDBS      = [];
dataDBS(:,1) = Cent_filt{m,1}(:,x); % x in mum
dataDBS(:,2) = Cent_filt{m,1}(:,y); % y in mum


if isempty(dataDBS)
   Cent_filt{m,6} = [];
else

% Run DBSCAN on each particle 

k   = 10;                                                   % minimum number of neighbors within Eps
Eps = 50;                                                   % minimum distance between points, nm

[class,type]=DBSCAN(dataDBS,k,Eps);                         % uses parameters specified at input
class2=transpose(class);                                    % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                                      % (core: 1, border: 0, outlier: -1)

coreBorder = [];
coreBorder = find(type2 >= 0);

subset          = [];
subset          = Cent_filt{m,1}(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);


% subsetP = [];
% 
% subsetP(:,1)    = dataDBS(coreBorder,1);
% subsetP(:,2)    = dataDBS(coreBorder,2);
% subsetP(:,3)    = class2(coreBorder);
% 
% figure('Position',[700 600 900 400])
% subplot(1,2,1)
% scatter(dataDBS(:,1),dataDBS(:,2),1);
% title('Raw Data')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on
% 
% subplot(1,2,2)
% scatter(subsetP(:,1),subsetP(:,2),1,mod(subsetP(:,3),10))
% title('identified Clusters')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on


% Select only the largest cluster(s)

if isempty(subset);
else

ClusterLength = [];

for i=1:max(subset(:,end)); % find the i-th cluster
    vx = find(subset(:,end)==i);
    ClusterLength(i,1) = length(vx);
end


% 6. Filtered data
% 7. Rg
% 8. Ecc

MaxClusterLength = find(ClusterLength(:,1) == max(ClusterLength));

vx = find(subset(:,end)==MaxClusterLength(1,1));

Cent_filt{m,6}  = subset(vx,1:end);

% Radius of Gyration equals the sum of the variances of x,y,z divided by
% the number of locs

Cent_filt{m,7}  = sqrt(sum(var(Cent_filt{m, 6}(:,x:y),1,1))); % Rg

% Eccentricity 
% covariance of x and y --> sqrt of min/max(Eigenvalues)

Cent_filt{m,8}  = sqrt(max(eig(cov(Cent_filt{m, 6}(:,x:y))))/min(eig(cov(Cent_filt{m, 6}(:,x:y))))); % Ecc


end
% for i=1:max(subset(:,end)); % find the i-th cluster
%     
%     vx = find(subset(:,end)==i);
%     
%     if length(vx)>70;
% 
%     Cent_filt{m,6}{i,1} = subset(vx,1:end);
%          
%     else    end
%     
% end


% Select only large cluster(s)

% for i=1:max(subset(:,end));
%     
%     vx = find(subset(:,end)==i);
%     
%     if length(vx)>70;
% 
%     Cent_filt{m,6}{i,1} = subset(vx,1:end);
%          
%     else    end
%     
% end

end

end

% figure('Position',[10 600 400 400],'name','# of Locs vs. GFP Intensity after DBSCAN');
% 
% for i=1:length(Cent_filt);
%     
%     for j = 1:length(Cent_filt{i,6});
%     
%     if isempty(Cent_filt{i,6}{j,1})
%     else
%     
%         scatter(length(Cent_filt{i,6}{j,1}),Cent_filt{i,2},5,'filled','r');hold on;
%         xlabel('Nbr of localizations');
%         ylabel('WF GFP integrated intensity');
%         box on;
%         axis square;
%     end
%     end
%     
% end

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)


%% Plot Rg and Ecc in subplots

% Set Threshold for Rg 

figure('Position',[500 300 1400 400])

subplot(1,3,1)
hist(cell2mat(Cent_filt(:,7)),30);
title('Radius of Gyration');
box on

subplot(1,3,2)
hist(cell2mat(Cent_filt(:,8)),30);
title('Eccentricity');
box on

subplot(1,3,3)
scatter(cell2mat(Cent_filt(:,7)),cell2mat(Cent_filt(:,8)),3);
title('Rg vs Ecc');
box on

%% Apply a filter for Rg and Ecc

select1 = find((cell2mat(Cent_filt(:,7))>10) & (cell2mat(Cent_filt(:,7))<200) & (cell2mat(Cent_filt(:,8))>2) & (cell2mat(Cent_filt(:,8))<5)); 

Cent_selected = {};
for i = 1:length(select1);  
Cent_selected{i,1} = Cent_filt{select1(i),1};
Cent_selected{i,2} = Cent_filt{select1(i),2};
Cent_selected{i,3} = Cent_filt{select1(i),3};
Cent_selected{i,4} = Cent_filt{select1(i),4};
Cent_selected{i,5} = Cent_filt{select1(i),5};
Cent_selected{i,6} = Cent_filt{select1(i),6};
Cent_selected{i,7} = Cent_filt{select1(i),7};
Cent_selected{i,8} = Cent_filt{select1(i),8};
end


NbrSubplots = round(sqrt(length(select1)))+1;
count = 1;
figure('Position',[100 100 1500 1000])

for i = 1:length(select1);
    
    
    if max(Cent_filt{select1(i),1}(:,4)) == 0 %max(Cent_filt{m,1}(:,1))/max(Cent_filt{m,1}(:,2))>1    
    
    x = 2; y = 3;

    else
    
    x = 1; y = 2;

    end
    
       
        subplot(NbrSubplots,NbrSubplots,count);
        scatter(Cent_filt{select1(i),6}(:,x),Cent_filt{select1(i),6}(:,y),10,'filled','black');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);
        count = count+1;
        box on
        axis off
        
   
end

%% Render the filtered images 

pxlsize = 10;

cd('Z:\Christian-Sieben\data_HTP\2016-10-12_Yeast_Wt_Kog1_GFP\analysis\rendered');

% Determine the box size form the largest particle

im_size = [];

for i=1:length(Cent_selected);
    
    
        if max(Cent_selected{i,1}(:,4)) == 0 %max(Cent_filt{m,1}(:,1))/max(Cent_filt{m,1}(:,2))>1    
    
        x = 2; y = 3;

        else
    
        x = 1; y = 2;
        
        end
    
im_size(i,1) = round((max(Cent_selected{i,6}(:,y))-min(Cent_selected{i,6}(:,y)))/pxlsize);
im_size(i,2) = round((max(Cent_selected{i,6}(:,x))-min(Cent_selected{i,6}(:,x)))/pxlsize);
        
        
        
end

count=1;

for i = 1:length(Cent_selected);
    
    
        if max(Cent_filt{select1(i),1}(:,4)) == 0 %max(Cent_filt{m,1}(:,1))/max(Cent_filt{m,1}(:,2))>1    
    
        x = 2; y = 3;

        else
    
        x = 1; y = 2;

        end
    
    
        heigth =round((max(Cent_selected{i,6}(:,x)) - min(Cent_selected{i,6}(:,x)))/pxlsize);
        width = round((max(Cent_selected{i,6}(:,y)) - min(Cent_selected{i,6}(:,y)))/pxlsize);
        
        rendered = hist3([Cent_selected{i,6}(:,y),Cent_selected{i,6}(:,x)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;        
 
name32rendered  = ['image_10nm_large_Particles_rendered' num2str(i),'_',num2str(j) '.tiff'];

% name16          = ['image_10nm_16bit_' num2str(i),'_',num2str(j) '.tiff'];    
% name32          = ['image_10nm_32bit_' num2str(i),'_',num2str(j) '.tiff'];

% I16 = uint16(empty);
% imwrite(I16, name16);


I32 = [];
% I32 = uint32(imgaussfilt(empty,1));
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


