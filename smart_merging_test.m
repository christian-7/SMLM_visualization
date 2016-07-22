%% Load Data

filename_peaks='humCent_aTubNB_Sas6_Pos1_1_MMStack_locResults_DC'; % filename of TS output file

filename_peaks2=[filename_peaks '.dat'];
peaks=dlmread(filename_peaks2,',',1,0);

file = fopen(filename_peaks2);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
frame = strmatch('frame',h);
photons = strmatch('intensity [photon]',h);
sigma = strmatch('sigma [nm]',h);

fprintf('\n -- Data Loaded --\n')

%% ROI

xmin=1.92*1e4;
xmax=2.0*1e4;

ymin=8.0*1e4;
ymax=8.06*1e4;


vx=find(peaks(:,x)>xmin & peaks(:,x)<xmax);
subset1=peaks(vx,1:end);
vy=find(subset1(:,y)>ymin & subset1(:,y)<ymax);
subset2=subset1(vy,1:end);

figure
% scatter(subset2(:,x),subset2(:,y),1,subset2(:,frame),'filled')
scatter(subset2(:,x),subset2(:,y),3,'filled');
xlabel('X (nm)');
ylabel('X (nm)');
box on;

length(subset2);

%% Create the input pos_list for the tracker

load('data_set2.mat')

pos_list(:,1)=all(:,1)+500;                   % in pxl
pos_list(:,2)=all(:,2)+500;                   % in pxl
pos_list(:,3)=all(:,5);             % photons
pos_list(:,4)=all(:,4);               % dt in seconds
pos_list=sortrows(pos_list,4);
%% 


pos_list(:,1)=subset2(:,x);                   % in pxl
pos_list(:,2)=subset2(:,y);                   % in pxl
pos_list(:,3)=subset2(:,photons);             % photons
pos_list(:,4)=subset2(:,frame);               % dt in seconds

%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

max_disp    = 15;           % in unit of data
min_pos     = 1;            % good - eliminate if fewer than good valid positions
gap         = 500;          % mem - number of time steps that a particle can be 'lost' and then recovered again
quiet       = 1;            % quiet - 1 = no text


param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res=trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

%% Plot All Tracks

% scatter(res(:,1),res(:,2),5,res(:,5));
% 
% for i=1:max(res(:,5))
%     
%     vx=find(res(:,5)==i);
%     
%     plot(res(vx,1),res(vx,2));hold on;
% 
% end

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
        
        track(j,6)= sqrt(((track_center(:,1)-track(j,1))^2)+((track_center(:,2)-track(j,2))^2)); % calculate distance from track_centers
        track(j,7)= track(j,4)-track(1,4); % time gap in frames
        
    end
    
    else end
    
    res(vx,6)=track(:,6); % distance
    res(vx,7)=track(:,7); % gap
  
end

% res1 = x
% res2 = y
% res3 = photons
% res4 = frame
% res5 = track ID
% res6 = distance
% res7 = gap time


%% Assign Probability to each localization per cluster

% /home/csiebenpc5/Documents

% load the fitting results

load('/Users/Christian/Documents/Arbeit/MatLab/centriole_sim/smart_merging/fit_gap_scatter.mat')


res(:,8)=zeros(1);      % total probability = probablity from distance * probablity from time

for i=1:length(res);  
   
    res(i,8)=fitresult(res(i,6))*gapfit(res(i,7));

end

% Show a histogram of probability and a scatter with color corresponding to
% probability

% figure ('Position',[100 500 800 300])
% 
% subplot(1,2,1)
% hist(res(:,8),20);
% % title('Total probability');
% xlabel('probability');
% ylabel('counts');
% 
% subplot(1,2,2)
% scatter(res(:,1),res(:,2),10,res(:,8),'filled')
% colorbar
% box on
% xlabel('x (nm)');
% ylabel('y (nm)');  

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

%% Merging/Weighting points

% Filter by probability

prob_tresh=0.003;
above_thresh=[];
below_thresh=[];

filter1=find(res(:,8)>prob_tresh);
filter2=find(res(:,8)<prob_tresh);

above_thresh=res(filter1,1:end);
below_thresh=res(filter2,1:end);

%% Merge Molecule above treshold

% res1 = x
% res2 = y
% res3 = photons
% res4 = frame
% res5 = track ID
% res6 = distance
% res7 = gap time
% res8 = probability to be a loc from the track


groupedx=[];
groupedy=[];
frame=[];
groupedframe=[];
groupedID=[];
groupedPhotons=[];
Photons=[];
above_thresh_merged=[];

for index=min(above_thresh(:,5)):max(above_thresh(:,5)); 
    
            vx=find(above_thresh(:,5)==index);
    
            clusterx=[];
            clustery=[];
            clusterxC=[];
            clusteryC=[];
                                                   
            clusterx=above_thresh(vx,1);
            clustery=above_thresh(vx,2);
            frame=above_thresh(vx,4);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            Photons=sum(above_thresh(vx,3));
            
            
            groupedx=vertcat(groupedx,clusterxC);
            groupedy=vertcat(groupedy,clusteryC);
            groupedframe=vertcat(groupedframe, round(mean(frame)));
            groupedID=vertcat(groupedID, index);
            groupedPhotons=vertcat(groupedPhotons, Photons);

end


above_thresh_merged(:,1)=groupedx(~isnan(groupedx));
above_thresh_merged(:,2)=groupedy(~isnan(groupedy));
above_thresh_merged(:,3)=groupedPhotons(~isnan(groupedx));
above_thresh_merged(:,4)=groupedframe(~isnan(groupedframe));;
above_thresh_merged(:,5)=groupedID(~isnan(groupedx));;
above_thresh_merged(:,6)=1;


fprintf('\n -- Localization above threshold merged --\n')

%% Join Merged Molecules and Locs below threshold

all_locs=zeros(length(above_thresh_merged)+length(below_thresh),4);

% x,y,photons,probablity

all_locs(:,1) = [above_thresh_merged(:,1); below_thresh(:,1)]; 
all_locs(:,2) = [above_thresh_merged(:,2); below_thresh(:,2)]; 
all_locs(:,3) = [above_thresh_merged(:,3); below_thresh(:,3)]; 
all_locs(:,4) = [above_thresh_merged(:,6); below_thresh(:,8)]; 


%% Determine the bin for each localization and wheight the pixel acc to the probability
% For the TH case

pxlsize=10;

heigth=round((max(all_locs(:,2))-min(all_locs(:,2)))/pxlsize);
width=round((max(all_locs(:,1))-min(all_locs(:,1)))/pxlsize);

% Generate the 2D Histogram for the thresholded data TH

imTHMerging = hist3([all_locs(:,1),all_locs(:,2)],[width heigth]); 

% Find the pixel for each localization

bin=[];

for i=1:length(all_locs); % for all molecules, find x and y pixel position
    
%     bin(i,1) = ceil(all_locs(i,1)/pxlsize)-ceil(min(all_locs(:,1))/pxlsize);
      bin(i,1) = round((all_locs(i,1)/pxlsize)-(min(all_locs(:,1))/pxlsize));
    
    if bin(i,1) == 0;
       bin(i,1) = 1;
    else end
    
    
%     bin(i,2) = round(all_locs(i,2)/pxlsize)-ceil(min(all_locs(:,2))/pxlsize);
      bin(i,2) = round((all_locs(i,2)/pxlsize)-(min(all_locs(:,2))/pxlsize));
    
     if bin(i,2) == 0;
       bin(i,2)=1;
    else end
    
    bin(i,3) = all_locs(i,3);%*all_locs(i,4); % Photons * Probability
    
end

% Calculate the number of photons in each pixel

newImage=[];

for i = 1:max(bin(:,1));
    
    for j = 1:max(bin(:,2));
    
    target=find(bin(:,1)==i & bin(:,2)==j);
    
    newImage(i,j)=sum(bin(target,3));
    
    end
    
end    


imTHMerging2 = times(imTHMerging, newImage);

%% Determine the bin for each localization and wheight the pixel acc to the probability
% For the Merged case

pxlsize=10;

heigth=round((max(merged(:,2))-min(merged(:,2)))/pxlsize);
width=round((max(merged(:,1))-min(merged(:,1)))/pxlsize);

% Generate the 2D Histogram for the smart merged data

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

imWMerging2 = times(imWMerging, newImage2);

%% Calculate 2D Histogram

% Apply Gaussian Filter to unmerged data

imUnmerged = hist3([pos_list(:,1),pos_list(:,2)],[width heigth]); 
G = fspecial('gaussian',[4 4],1.5); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imUnmerged_G = imfilter(imUnmerged,G,'same');

% Apply Gaussian filter to TH 2D histogram (TH before)

G = fspecial('gaussian',[4 4],1.5); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imTHMerging_G = imfilter(imTHMerging,G,'same');

% Apply Gaussian Filter to weighted 2D histogram (TH after)

imTHMerging2_G = imfilter(imTHMerging2,G,'same');

% Apply Gaussian Filter to the smart merged  2D histogram (SM before)

G = fspecial('gaussian',[4 4],1.5); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imWMerging_G = imfilter(imWMerging,G,'same');

% Apply Gaussian Filter to the smart merged  2D histogram (SM after)

G = fspecial('gaussian',[4 4],1.5); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imWMerging2_G = imfilter(imWMerging2,G,'same');

% Show 2D 

figure%('Position',[100 400 1000 600])
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[100 400 900 600]);

subplot(2,3,1)
imagesc(imUnmerged_G);
title('unmerged')
colormap('hot')

subplot(2,3,2)
imagesc(imTHMerging_G)
title('merged with TH')
colormap('hot');

subplot(2,3,3)
imagesc(imTHMerging2_G);
title('merged with TH + weighting')
colormap('hot');% colorbar;

subplot(2,3,4)
imagesc(imWMerging_G);
title('weighted merging')
colormap('hot'); %colorbar;

subplot(2,3,5)
imagesc(imWMerging2_G);
title('weighted merging + weighting')
colormap('hot'); %colorbar;

%% Show QuadTree rendering without photon count weighting
% 
% [ I1 ] = QuadTree( subset2(:,1:2), 15,255);
% [ I2 ] = QuadTree( all_locs(:,1:2), 15,255);
% [ I3 ] = QuadTree( merged(:,1:2), 10,255);
% 
% figure('Position', [100 100 900 300]);
% 
% subplot(1,3,1)
% imagesc(I1, [1 30]);
% title('Unmerged');
% colormap('hot');
% 
% subplot(1,3,2);
% imagesc(I2,[1 50]);
% title('After TH');
% colormap('hot');
% 
% subplot(1,3,3);
% imagesc(I3,[1 50]);
% title('Merged');
% colormap('hot');


%% Write
imwrite(imWMerging2_G,'imWMerging2_G.tiff');
imwrite(imTHMerging2_G,'imTHMerging2_G.tiff');

%% Write 32-bit Files
tic
% % Write 32-bit file

% cd('Z:\Christian-Sieben\data_HTP\2016-02-18_Centriole_Nanobody\locResults\rendered_matlab');
I32=[];
I32=uint32(imWMerging2_G);

outputFileName1 = ['Smart_Merged_Gaussian' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

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

%%%

I32=[];
I32=uint32(imWMerging2);

outputFileName1 = ['Smart_Merged_' num2str(pxlsize) 'nm_pxl_32bit.tiff'];

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

fprintf(' -- Images saved (%f sec) -- \n',toc)




