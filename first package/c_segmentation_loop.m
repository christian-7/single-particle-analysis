%% Loop for centriole segmentation

% base = 'FOV1_Sas6_1000mW_10ms_A647_';

for i = [8,9,10];
    
    base = ['FOV' num2str(i) '_Sas6_1000mW_10ms_A647_1'];

% Set file and path names

% WFpath = ['Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\humanCent_Sas6_A647_WF' num2str(i)];
WFpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647';
% WF_name = [base 'WF' num2str(i) '_MMStack_Pos0.ome.tif'];
WF_name = ['SAS6_FOV' num2str(i) '.tif'];

Locpath = ['Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\' base];
locName = [base '_MMStack_Pos0_locResults_DC.dat'];

savepath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles';
% savename = [base num2str(i), '_segmented.mat'];
savename = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(i), '_segmented.mat'];

% Start segmentation

c_image_segmentation(WFpath,WF_name,Locpath,locName,savepath,savename);

fprintf('\n -- %f finished --\n',i)

end

%% Loop for smart merging of the segmented centrioles

base = 'humanCent_Sas6_A647_';

for m = 1:10;

% Locpath = ['Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\locResults_Sas6\humanCent_Sas6_A647_' num2str(m)];
% locName = [base num2str(m),'_MMStack_Pos0_locResults_DC.dat'];

Locpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\FOV1_Sas6_1000mW_10ms_A647_1';
locName = 'FOV1_Sas6_1000mW_10ms_A647_1_MMStack_Pos0_locResults_DC.dat';

Segpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles';
segName = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(m), '_segmented.mat'];

c_smart_merging(Locpath,locName,Segpath,segName)

end

%% Loop for image generation

base = 'humanCent_Sas6_A647_';
pxlsize = 10;

for m = 1:10;

% Locpath = ['Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\locResults_Sas6\humanCent_Sas6_A647_' num2str(m)];
% locName = [base num2str(m),'_MMStack_Pos0_locResults_DC.dat'];

Locpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\FOV1_Sas6_1000mW_10ms_A647_1';
locName = 'FOV1_Sas6_1000mW_10ms_A647_1_MMStack_Pos0_locResults_DC.dat';

Segpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles';
% segName =  [base num2str(m), '_segmented.mat'];
segName = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(m), '_segmented.mat'];

% savename = [base num2str(m), '_smart_merged.mat'];    
savename = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(m), '_smart_merged.mat'];  % name for the concatenated list of smart merged locs

impath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles\2D_hist_after_smart_merging';  % path for the 32 bit generated 2D Hostrogram
% outputFileName1 = [base num2str(m),'_smart_merged_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];                                      % name for the 2D Histogram
outputFileName1 = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(m) '_smart_merged_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

c_render_image(Locpath,locName,Segpath,segName,savename,impath,outputFileName1);

end

%% Loop for generating an image without merging

base = 'humanCent_aTub_NB_A647_';
pxlsize = 10;

for m = [10];

% Locpath = ['Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\locResults_aTub\humanCent_aTub_NB_A647_' num2str(m)];
% locName = [base num2str(m),'_MMStack_Pos0_locResults_DC.dat'];

Locpath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\FOV10_Sas6_1000mW_10ms_A647_1';
locName = 'FOV10_Sas6_1000mW_10ms_A647_1_MMStack_Pos0_locResults_DC.dat';

impath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\images_without_merging';  % path for the 32 bit generated 2D Hostrogram
% outputFileName1 = [base num2str(m),'_', num2str(pxlsize) '_nm_pxl_32bit.tiff'];                                      % name for the 2D Histogram
outputFileName1 = ['humanCent_Sas6_1000mW_10ms_A647_' num2str(m) '_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

c_render_image_wo_merging(Locpath,locName,pxlsize,impath,outputFileName1) 

fprintf('\n -- %f finished --\n',m)

end




