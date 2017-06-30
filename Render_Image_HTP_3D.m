clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

savename        = 'Tom20_A647_10ms_5';
savepath        = 'Z:\Christian-Sieben\data_HTP\2017-05-22_3D_Test_Mito\locResults\Tom20_A647_10ms_5';

locpath_Ch1     = 'Z:\Christian-Sieben\data_HTP\2017-05-22_3D_Test_Mito\locResults\Tom20_A647_10ms_5';
locname_Ch1     = 'Tom20_A647_10ms_5_MMStack_1_Localizations_Z_DC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(locpath_Ch1);
locs_Ch1        = dlmread([locname_Ch1 '.csv'],',',1,0);

cd(locpath_Ch1);
file = fopen([locname_Ch1 '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol                = strmatch('x [nm]',h);
yCol                = strmatch('y [nm]',h);
zCol                = strmatch('z [nm]',h);
framesCol           = strmatch('frame',h);
LLCol               = strmatch('loglikelihood',h);
photonsCol          = strmatch('intensity [photon]',h);
sigmaXCol           = strmatch('sigma_x [nm]',h);
sigmaYCol           = strmatch('sigma_y [nm]',h);
uncertaintyCol      = strmatch('uncertainty [nm]',h);

xCol                = strmatch('"x [nm]"',h);
yCol                = strmatch('"y [nm]"',h);
zCol                = strmatch('"z [nm]"',h);
framesCol           = strmatch('"frame"',h);
LLCol               = strmatch('"loglikelihood"',h);
photonsCol          = strmatch('"intensity [photon]"',h);
sigmaXCol           = strmatch('"sigma_x [nm]"',h);
sigmaYCol           = strmatch('"sigma_y [nm]"',h);
uncertaintyCol      = strmatch('"uncertainty [nm]"',h);

fprintf('\n -- Data Loaded --\n')

%% Show filtering parameters

close all;

filter2 = find(locs_Ch1(:,zCol)==0);
locs_Ch1(filter2,:) = [];

figure('Position',[100 100 800 500]); 
subplot(2,3,1);
hist(locs_Ch1(locs_Ch1(:,photonsCol)<2e4,photonsCol),50);
title('Photons');
subplot(2,3,2);
hist(locs_Ch1(:,sigmaXCol),50);
title('Sigma X');
subplot(2,3,3);
hist(locs_Ch1(locs_Ch1(:,LLCol)<1500,LLCol),50);
title('LL');
subplot(2,3,4);
hist(locs_Ch1(locs_Ch1(:,uncertaintyCol)<100,uncertaintyCol),50);
title('Uncertainty');
subplot(2,3,5);
hist(locs_Ch1(:,sigmaYCol),50);
title('Sigma Y');
subplot(2,3,6);
hist(locs_Ch1(:,zCol),50);
title('Z Range');

%% Choose and apply the filter parameters

close all

% Set Filter parameters for Channel 1

minFrame            = 1000;
MinPhotons          = 200;
MinZ                = -500; 
MaxZ                = 500;
Maxuncertainty      = 30;
logLikelihood       = 1000;

filter   = find(locs_Ch1(:,zCol) < MaxZ & locs_Ch1(:,zCol) > MinZ & locs_Ch1(:,photonsCol) > MinPhotons & locs_Ch1(:,uncertaintyCol) < Maxuncertainty & locs_Ch1(:,framesCol) > minFrame & locs_Ch1(:,LLCol) > logLikelihood);
locs_Filtered_Ch1 = locs_Ch1(filter,1:end);

fprintf('\n -- Ch1 Filtered %f of the locs are left --\n', ((length(locs_Filtered_Ch1)/length(locs_Ch1))));

%% Select an ROI

% Plot the data
% Show 2D histogram 
% Select Au Fuducial using rectangular selection

close all
% 
% figure('Position',[100 400 500 500])
% scatter(locs(:,xcol),locs(:,ycol),1);

pxlsize = 50; 

heigth  = round((max(locs_Filtered_Ch1(:,yCol))-min(locs_Filtered_Ch1(:,yCol)))/pxlsize);
width   = round((max(locs_Filtered_Ch1(:,xCol))-min(locs_Filtered_Ch1(:,xCol)))/pxlsize);

figure('Position',[650 400 500 500])
im=hist3([locs_Filtered_Ch1(:,xCol),locs_Filtered_Ch1(:,yCol)],[width heigth]); % heigth x width
imagesc(imrotate(im,90),[0 5]);
% imagesc(im,[0 200]);
colormap('hot');
% colorbar
rect = getrect; % rect = [xmin ymin width height];

close all

fprintf('\n -- ROI selected --\n')


%% Select ROI 
close all

xmin = min(locs_Filtered_Ch1(:,xCol))+ rect(:,1)*pxlsize;
ymin = max(locs_Filtered_Ch1(:,yCol)) - rect(:,2)*pxlsize - (rect(:,4)*pxlsize) ;
xmax = xmin + (rect(:,3)* pxlsize);
ymax = ymin + rect(:,4) * pxlsize;

% xmin=min(locs(:,3));
% xmax=max(locs(:,3)); 
% ymin=min(locs(:,4));
% ymax=max(locs(:,4));

vx      = find(locs_Filtered_Ch1(:,xCol)>xmin & locs_Filtered_Ch1(:,xCol)<xmax & locs_Filtered_Ch1(:,yCol)>ymin & locs_Filtered_Ch1(:,yCol)<ymax);
subset1 = locs_Filtered_Ch1(vx,1:end);

figure('Position',[1200 400 500 500])
scatter(subset1(:,xCol),subset1(:,yCol),2,subset1(:,zCol))
colormap('jet');
colorbar

figure('Position',[1200 400 500 500])
scatter3(subset1(:,xCol),subset1(:,yCol),subset1(:,zCol),2,'filled');hold on;
colorbar

locs_Filtered_Ch1 = subset1;

fprintf('\n -- Plotted selected ROI  --\n')

%% Export ROI for VISP

% VISP format x,y,z,Photons, Frame

cd(locpath_Ch1);

locs_subset = [];

locs_subset (:,1) = locs_Filtered_Ch1(:,xCol);
locs_subset (:,2) = locs_Filtered_Ch1(:,yCol);
locs_subset (:,3) = locs_Filtered_Ch1(:,zCol);
locs_subset (:,4) = locs_Filtered_Ch1(:,photonsCol);
locs_subset (:,5) = locs_Filtered_Ch1(:,framesCol);

name_for_VISP = [savename, '_ROI2_forVISP.txt'];

% fid = fopen(name_for_LAMA,'wt');
% fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_VISP, locs_subset,'delimiter', '\t');
% fclose(fid);


fprintf('\n -- Data saved in VISP format --\n')





%% Determine slice for each localization

SliceThickness      = 100; % nm
slice               = zeros(length(locs_Filtered_Ch1),1);
number_of_slices    = round((max(locs_Filtered_Ch1(:,zCol)) - min(locs_Filtered_Ch1(:,zCol)))/SliceThickness);

for i=1:max(number_of_slices);
    
    selectedSlice = find(locs_Filtered_Ch1(:,zCol) > min(locs_Filtered_Ch1(:,zCol))+((i-1)*SliceThickness) & locs_Filtered_Ch1(:,zCol) < min(locs_Filtered_Ch1(:,zCol))+(i*SliceThickness));
    
    slice(selectedSlice,1) = i;
         
end

number_of_slices

fprintf('\n -- Determined slice for each localization --\n');

%% Build Box around each slice

xmin = 0;
ymin = 0;
xmax = max(locs_Filtered_Ch1(:,xCol));
ymax = max(locs_Filtered_Ch1(:,yCol));

interval = 10;

lowerLine = [];
lowerLine(:,1) = xmin:((max(xmax))/interval):max(xmax); 
lowerLine(1:end,2) = ymin; 

upperLine = [];
upperLine(:,1) = xmin:(max(xmax)/interval):max(xmax); 
upperLine(1:end,2) = max(ymax)*1.01; 

leftLine = [];
leftLine(:,2) = ymin:((max(ymax))/interval):max(ymax); 
leftLine(1:end,1) = xmin; 

rightLine = [];
rightLine(:,2) = ymin:(max(ymax)/interval):max(ymax); 
rightLine(1:end,1) = max(xmax)*1.01; 

addedBox = [lowerLine; upperLine; leftLine; rightLine];

figure
scatter(addedBox(:,1), addedBox(:,2));


%% Calculate 2D histogram for each slice

pxlsize = 20;

for i = 1:max(slice(:,1));
    
selectedSlice = locs_Filtered_Ch1(slice(:,1)==i,1:end);
    
selectedSlice = [selectedSlice(:,xCol:yCol);addedBox];
    
heigth  = round((max(selectedSlice(:,2))-min(selectedSlice(:,2)))/pxlsize);
width   = round((max(selectedSlice(:,1))-min(selectedSlice(:,1)))/pxlsize);

% heigth  = round((max(addedBox(:,yCol))-min(addedBox(:,yCol)))/pxlsize);
% width   = round((max(addedBox(:,xCol))-min(addedBox(:,xCol)))/pxlsize);

im = hist3([selectedSlice(:,2),selectedSlice(:,1)],[heigth width]);

% % imG = imgaussfilt(im, pxlsize/10);

cd(savepath)

name = [savename '_rendered_' num2str(pxlsize) 'nm_per_pxl_slice_' num2str(i) '.tiff'];  

I32=[];
I32=uint32(im);

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

end

fprintf('\n -- Saved rendered slices --\n');

%% Save for VISP

% VISP format x,y,z, Photons, Frame

locs_subset = [];

locs_subset (:,1) = locs_Ch1(filter,xCol);
locs_subset (:,2) = locs_Ch1(filter,yCol);
locs_subset (:,3) = locs_Ch1(filter,zCol);
locs_subset (:,4) = locs_Ch1(filter,photonsCol);
locs_subset (:,5) = locs_Ch1(filter,framesCol);

name_for_VISP = [locname_Ch1, '_forVISP.txt'];

% fid = fopen(name_for_LAMA,'wt');
% fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

cd(savepath)
dlmwrite(name_for_VISP, locs_subset,'delimiter', '\t');
% fclose(fid);


fprintf('\n -- Data saved in VISP format --\n')


