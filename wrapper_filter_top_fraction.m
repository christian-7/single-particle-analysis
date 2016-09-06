%% Wrapper

% all_particles = [];
% count=1;

for i=1:length(FOV10_Sas6_1000mW_10ms_A647_1_extractedParticles);
    
    for j = 1:length(FOV10_Sas6_1000mW_10ms_A647_1_extractedParticles{i,4});
        
        if isempty(FOV10_Sas6_1000mW_10ms_A647_1_extractedParticles{i,4}{j,1})
        else
    all_particles{count,1} = FOV10_Sas6_1000mW_10ms_A647_1_extractedParticles{i,4}{j,1};
    count = count+1;
        end
        end
end


%% Select only the top 10 %

locs = [];

for i=1:length(all_particles);
    
locs(i,1) = length(all_particles{i,1});

end

% select first 10%

locs = sort(locs,'ascend');

threshold = locs(round((length(locs/100)*0.9)));


top_particles = [];
count=1;

for i=1:length(all_particles);
    
    if length(all_particles{i,1})>threshold;
    
    top_particles{count,1} = all_particles{i,1};
    count = count+1;
    else
    end
end

%% Plot Top_particles

nbrSubplots = round(sqrt(length(top_particles)));

figure('Position',[100 100 1500 1000]);

for i = 1:length(top_particles);
   
    subplot(nbrSubplots+1,nbrSubplots,i)
    scatter(top_particles{i,1}(:,1),top_particles{i,1}(:,2),1,'black');
    title(num2str(i));
    axis off  
    
end

%% Plot Top_particles

pxlsize = 5;
cd('Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles\extracted_DBSCAN_filtered\rendered');

% Determine the box size form the largest particle

im_size = [];



for i = 1:length(top_particles);


im_size(i,1) = round((max(top_particles{i,1}(:,2))-min(top_particles{i,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(top_particles{i,1}(:,1))-min(top_particles{i,1}(:,1)))/pxlsize);
        

end


count=1;

for i = [2,6,8,10,14,19,26,37,36,48,58] % 1:length(top_particles);
      
            
        heigth =round((max(top_particles{i,1}(:,1)) - min(top_particles{i,1}(:,1)))/pxlsize);
        width = round((max(top_particles{i,1}(:,2)) - min(top_particles{i,1}(:,2)))/pxlsize);
        
        rendered = hist3([top_particles{i,1}(:,2),top_particles{i,1}(:,1)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;
 

    
name = ['image_5nm_filtered' num2str(count) '.tiff'];

I32=[];
I32=uint32(empty);

t = Tiff(name,'w');

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



