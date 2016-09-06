

pxlsize = 5;
cd('Y:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\FOV1_Sas6_1000mW_10ms_A647_1\rendered');

% Determine the box size form the largest particle

im_size = [];

for i=1:length(Cent);
    
    for j = 1:length(Cent{i,4});
        
        if isempty(Cent{i,4}{j,1})
        else

im_size(i,1) = round((max(Cent{i,4}{j,1}(:,2))-min(Cent{i,4}{j,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(Cent{i,4}{j,1}(:,1))-min(Cent{i,4}{j,1}(:,1)))/pxlsize);
        
        end
    end
end

count=1;

for i = 1:length(Cent);
    
    for j = 1:length(Cent{i,4});
    
            if isempty(Cent{i,4}{j,1})
                
            else
            
        heigth =round((max(Cent{i,4}{j,1}(:,1)) - min(Cent{i,4}{j,1}(:,1)))/pxlsize);
        width = round((max(Cent{i,4}{j,1}(:,2)) - min(Cent{i,4}{j,1}(:,2)))/pxlsize);
        
        rendered = hist3([Cent{i,4}{j,1}(:,2),Cent{i,4}{j,1}(:,1)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;
 

    
name = ['image_10nm_' num2str(count) '.tiff'];

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
    end
end
        