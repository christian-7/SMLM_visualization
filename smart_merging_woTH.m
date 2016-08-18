clc, clear, close all

%% Load Data

cd('.\simulated_test_data'); 
peaks = dlmread('suresim_1localizations.txt',' ',1,0);
xCol = 1;
yCol = 2;
frameCol = 4;
photonsCol = 5;

% cd('.\simulated_test_data'); 
% peaks = sim_line;
% xCol = 1;
% yCol = 2;
% frameCol = 4;
% photonsCol = 3;

% cd('.\test_data'); 
% load('real_Cent_large.mat')
% peaks = subset2;
% 
% xCol = 1;%2
% yCol = 2;%3
% frameCol = 4;%5
% photonsCol = 6;%6

suffix = 'sureSim_1';


% filename_peaks2=[filename_peaks '.dat'];
% peaks=dlmread(filename_peaks2,',',1,0);

% file = fopen(filename_peaks2);
% line = fgetl(file);
% h = regexp( line, ',', 'split' );
% 
% x = strmatch('x [nm]',h);
% y = strmatch('y [nm]',h);
% frame = strmatch('frame',h);
% photons = strmatch('intensity [photon]',h);
% sigma = strmatch('sigma [nm]',h);
cd('..\')

fprintf('\n -- Data Loaded --\n')

%% ROI
% 
% xmin=1.92*1e4;
% xmax=2.0*1e4;
% 
% ymin=8.0*1e4;
% ymax=8.06*1e4;
% 
% 
% vx=find(peaks(:,x)>xmin & peaks(:,x)<xmax);
% subset1=peaks(vx,1:end);
% vy=find(subset1(:,y)>ymin & subset1(:,y)<ymax);
% subset2=subset1(vy,1:end);
% 
% figure
% % scatter(subset2(:,x),subset2(:,y),1,subset2(:,frame),'filled')
% scatter(subset2(:,x),subset2(:,y),3,'filled');
% xlabel('X (nm)');
% ylabel('X (nm)');
% box on;
% 
% length(subset2);

%% Create input for the tracker

peaks=sortrows(peaks,4);

pos_list(:,1)=peaks(:,xCol);                   % in pxl
pos_list(:,2)=peaks(:,yCol);                   % in pxl
pos_list(:,3)=peaks(:,photonsCol);             % photons
pos_list(:,4)=peaks(:,frameCol);               % dt in frames


%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)
tic

max_disp    = 15;           % in unit of data, 15
min_pos     = 1;            % good - eliminate if fewer than good valid positions
gap         = 50;          % mem - number of time steps that a particle can be 'lost' and then recovered again, 100
quiet       = 1;            % quiet - 1 = no text


param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res=trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done after %f sec--\n', toc)

% Output: 
% 
% x, y, photons, frames, ID

%% Direkt merging 

groupedx=[];
groupedy=[];
frame=[];
direct_merging=[];

for index=1:max(res(:,5));              % find the ID
    
            vx=find(res(:,5)==index);
    
            clusterx=[];
            clustery=[];
            clusterxC=[];
            clusteryC=[];
                                                   
            clusterx=res(vx,1);
            clustery=res(vx,2);
            frame=res(vx,4);
            
            clusterxC = sum(clusterx)/length(clusterx);
            clusteryC = sum(clustery)/length(clustery);         
            
            groupedx = vertcat(groupedx,clusterxC);
            groupedy = vertcat(groupedy,clusteryC);

end

direct_merging(:,1)=groupedx;
direct_merging(:,2)=groupedy;


fprintf('\n -- Direkt merging done --\n'); 

%% Calculate the distance and the gap time from origin of track

res(:,6)=zeros(1); % distance r from the center of the track
res(:,7)=zeros(1); % gap time between 1st and n-th localization

for i=1:max(res(:,5));          % for all tracks
    
    vx=find(res(:,5)==i);       % find the i-th track 

    track=res(vx,1:5);          % single track
    
    track_center(:,1) = sum(track(:,1))/length(track(:,1)); % center of the i-th track, x
    track_center(:,2) = sum(track(:,2))/length(track(:,2)); % center of the i-th track, y
    
    track(:,6)=zeros(1);
    track(:,7)=zeros(1);
    
    if length(track(:,1))>1;
    
    for j=1:length(track(:,1));
        
        track(j,6)= sqrt(((track_center(:,1)-track(j,1))^2)+((track_center(:,2)-track(j,2))^2));    % calculate distance from track_centers
        track(j,7)= track(j,4)-track(1,4);                                                          % time gap in frames
        
    end
    
    else end
    
    res(vx,6)=track(:,6); % distance
    res(vx,7)=track(:,7); % gap
  
end

% Updated variable res

% res1 = x
% res2 = y
% res3 = photons
% res4 = frame
% res5 = track ID
% res6 = distance
% res7 = gap time


%% Assign Probability to each localization per cluster
% load the fitting results

load('K:\Christian\GitHub\SMLM_vis\exp_dist\dist_gap_fit_exp.mat');

res(:,8)=zeros(1);      % total probability = probablity from distance * probablity from time

for i=1:length(res);  
   
    res(i,8)=distfit(res(i,6))*gap_fit_exp(res(i,7));

end

fprintf('\n -- Probability assigned --\n')

% Show a histogram of probability and a scatter with color corresponding to
% probability



figure ('Position',[100 500 900 300])

subplot(1,2,1)
hist(res(:,8)/max(res(:,8)),10);
title('Total probability');
xlabel('probability');
ylabel('counts');
axis square

subplot(1,2,2)
scatter(res(:,1),res(:,2),5,res(:,8)/max(res(:,8)),'filled')
% axis([0 1000 -50 100])
colorbar
box on
xlabel('x (nm)');
ylabel('y (nm)');  

%% Weighted Merging according to the probability

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
    

for i=1:max(res(:,5));          % for all tracks
    
    vx=find(res(:,5)==i);       % find the i-th track 

    track=res(vx,1:end);        % single track
    
    ProbFrac = (track(:,8)/(sum(track(:,8))))*100;  % calculate probability fraction (frac of 100)
    
    for j=1:length(vx);
       
    x(j) = track(j,1)*ProbFrac(j);  % x coordinate
    y(j) = track(j,2)*ProbFrac(j);  % y coordinate
    Averagex=sum(x)/100;            % calculate the average
    Averagey=sum(y)/100;            % calculate the average
    
    end
    
            groupedx=vertcat(groupedx,Averagex);
            groupedy=vertcat(groupedy,Averagey);
            groupedframe=vertcat(groupedframe, round(mean(track(:,4))));
            groupedID=vertcat(groupedID, i);
            groupedPhotons=vertcat(groupedPhotons, sum(track(:,3)));
            
    x = [];
    y = [];
    Averagex = [];
    Averagey = [];
    
    
end

% x,y,photons,frame, ID

merged(:,1)=groupedx;
merged(:,2)=groupedy;
merged(:,3)=groupedPhotons;
merged(:,4)=groupedframe;
merged(:,5)=groupedID;

%% Determine the bin for each localization and wheight the pixel acc to the probability
% For the Samrt Merged case

pxlsize = 5;

heigth = round((max(merged(:,2))-min(merged(:,2)))/pxlsize);
width =  round((max(merged(:,1))-min(merged(:,1)))/pxlsize);

% Generate the 2D Histogram for the smart merged data, before photon weighting

imWMerging = hist3([merged(:,1),merged(:,2)],[width heigth]); 

% Find the pixel for each localization

binMerged=[];

for i=1:length(merged); % for all molecules, find x and y pixel
    
%   binMerged(i,1) = ceil(merged(i,1)/pxlsize)-ceil(min(merged(:,1))/pxlsize);
    binMerged(i,1) = round((merged(i,1)/pxlsize)-(min(merged(:,1))/pxlsize));
    
    if binMerged(i,1) == 0;
       binMerged(i,1)=1;
    else end
    
    
%     binMerged(i,2) = ceil(merged(i,2)/pxlsize)-ceil(min(merged(:,2))/pxlsize);
      binMerged(i,2) = round((merged(i,2)/pxlsize)-(min(merged(:,2))/pxlsize));
      
     if binMerged(i,2) == 0;
       binMerged(i,2)=1;
    else end
    
    binMerged(i,3) = merged(i,3); % *all_locs(i,4); % Photons * Probability
    
end

% Calculate the number of photons in each pixel

newImage2 = [];

for i = 1:max(binMerged(:,1));
    
    for j = 1:max(binMerged(:,2));
    
    target=find(binMerged(:,1)==i & binMerged(:,2)==j);
    
    newImage2(i,j)=sum(binMerged(target,3));
    
    end
    
end    

imWMerging2 = times(imWMerging, newImage2); % smart merged after photon weighting

%% Calculate 2D Histogram

Gfilter = 1.5;

% Apply Gaussian Filter to unmerged data

heigth = round((max(pos_list(:,2))-min(pos_list(:,2)))/pxlsize);
width  = round((max(pos_list(:,1))-min(pos_list(:,1)))/pxlsize);

imUnmerged = hist3([pos_list(:,1),pos_list(:,2)],[width heigth]); 
imUnmerged_G = imgaussfilt(imUnmerged,Gfilter);

% Apply Gaussian Filter to the smart merged  2D histogram (SM after)
imWMerging2_G = imgaussfilt(imWMerging2,Gfilter);

% Apply Gaussian Filter to the direct merged 2D histogram 
heigth = round((max(direct_merging(:,2))-min(direct_merging(:,2)))/pxlsize);
width  = round((max(direct_merging(:,1))-min(direct_merging(:,1)))/pxlsize);

dmerged = hist3([direct_merging(:,1),direct_merging(:,2)],[width heigth]); 
dmerged_G = imgaussfilt(dmerged,Gfilter);

% Show 2D 

figure
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[100 400 1500 150]);

subplot(1,3,1)
imagesc(imrotate(imUnmerged_G,90));
title('unmerged')
colormap('hot');
axis off

subplot(1,3,2)
imagesc(imrotate(dmerged_G,90));
title('direct merged')
colormap('hot'); %colorbar;
axis off

subplot(1,3,3)
imagesc(imrotate(imWMerging2_G,90));
title('smart merging + photon weighting')
colormap('hot'); %colorbar;
axis off

%% Write 8-bit tiff files

% suffix = 'real_MT_large_FOV_gap50';

imwrite(imUnmerged_G,[suffix '_unmerged' num2str(pxlsize) 'nm_pxl_8bit.tiff']);
imwrite(imWMerging2_G,[suffix '_smart_merged_' num2str(pxlsize) 'nm_pxl_8bit.tiff']);
imwrite(dmerged_G,[suffix '_direct_merging' num2str(pxlsize) 'nm_pxl_8bit.tiff']);



%% Write 32-bit File

% suffix = 'real_MT_2_gap_50';

I32=[];
I32=uint32(imUnmerged_G);

outputFileName1 = [suffix '_unmerged_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

t = Tiff(outputFileName1,'w');

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
t.close();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I32=[];
I32=uint32(imWMerging2_G);

outputFileName1 = [suffix '_smart_merged_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

t = Tiff(outputFileName1,'w');

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
t.close();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I32=[];
I32=uint32(dmerged_G);

outputFileName1 = [suffix '_direct_merging_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

t = Tiff(outputFileName1,'w');

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
t.close();

fprintf(' -- Images saved -- \n')




