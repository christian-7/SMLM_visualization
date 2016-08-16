%% Plot Data from SureSim

filename_peaks='suresim_1Localizations'; % filename of TS output file
filename_peaks2=[filename_peaks '.txt'];
peaks=dlmread(filename_peaks2,' ',1,0);

file = fopen(filename_peaks2);
line = fgetl(file);
h = regexp( line, ' ', 'split' );

xCol = strmatch('Pos_x',h);
yCol = strmatch('Pos_y',h);
frameCol = strmatch('Frame',h);
photonsCol = strmatch('Intensity',h);

fprintf('\n -- Data Loaded --\n')

% Plot the raw data

pxlsize=10;
sigma = 0.7; 

heigth=round((max(peaks(:,yCol))-min(peaks(:,yCol)))/pxlsize);
width=round((max(peaks(:,xCol))-min(peaks(:,xCol)))/pxlsize);

rendered_image = hist3([peaks(:,xCol),peaks(:,yCol)],[width heigth]); 
gauss_filt = imgaussfilt(rendered_image,sigma);
I32=uint32(gauss_filt);

figure('Position',[500 500 500 300])
subplot(1,2,1)
scatter(peaks(:,xCol),peaks(:,yCol),1,'black');
axis([min(peaks(:,xCol)) max(peaks(:,xCol)) min(peaks(:,yCol)) max(peaks(:,yCol))])
axis square
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(1,2,2)
imshow(imrotate(I32,90),[0 10]);
colormap hot
title('Blurred 2D Histogram');
axis square
axis off

%% Plot Data from SimpleSimulator

load('/Users/Christian/GitHub/SMLM_visualization/simulated_test_data/simulated_MT_3D_radius_30nm.mat')

peaks = sim_line;

xCol = 1;
yCol = 2;
frameCol = 4;
photonsCol = 3;

% Plot the raw data

pxlsize=10;
sigma = 0.7; 

heigth=round((max(peaks(:,yCol))-min(peaks(:,yCol)))/pxlsize);
width=round((max(peaks(:,xCol))-min(peaks(:,xCol)))/pxlsize);

rendered_image = hist3([peaks(:,xCol),peaks(:,yCol)],[width heigth]); 
gauss_filt = imgaussfilt(rendered_image,sigma);
I32=uint32(gauss_filt);

figure('Position',[500 500 500 300])
subplot(1,2,1)
scatter(peaks(:,xCol),peaks(:,yCol),1,'black');
axis([min(peaks(:,xCol)) max(peaks(:,xCol)) min(peaks(:,yCol)) max(peaks(:,yCol))])
axis square
axis([0 1e3 -50 100])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(1,2,2)
imshow(imrotate(I32,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

%% Plot Centriole Data from measurement

cent1 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_Cent_1.mat')
cent2 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_Cent_2.mat')
cent3 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_Cent_large.mat')

% x,y,photons, uncertainty, frame

peaks1 = cent1.subset2;
peaks2 = cent2.subset2;
peaks3 = cent3.subset2;

xCol = 1;
yCol = 2;
frameCol = 5;
photonsCol = 3;

% Plot the raw data for Cent1

pxlsize=10;
sigma = 0.7; 

heigth1=round((max(peaks1(:,yCol))-min(peaks1(:,yCol)))/pxlsize);
width1=round((max(peaks1(:,xCol))-min(peaks1(:,xCol)))/pxlsize);

rend1 = hist3([peaks1(:,xCol),peaks1(:,yCol)],[width1 heigth1]); 
gauss_filt_1 = imgaussfilt(rend1,sigma);
I32_1=uint32(gauss_filt_1);

% Plot the raw data for Cent2

heigth2=round((max(peaks2(:,yCol))-min(peaks2(:,yCol)))/pxlsize);
width2=round((max(peaks2(:,xCol))-min(peaks2(:,xCol)))/pxlsize);

rend2 = hist3([peaks2(:,xCol),peaks2(:,yCol)],[width2 heigth2]); 
gauss_filt_2 = imgaussfilt(rend2,sigma);
I32_2=uint32(gauss_filt_2);

% Plot the raw data for Cent3

heigth3=round((max(peaks3(:,yCol))-min(peaks3(:,yCol)))/pxlsize);
width3=round((max(peaks3(:,xCol))-min(peaks3(:,xCol)))/pxlsize);

rend3 = hist3([peaks3(:,xCol),peaks3(:,yCol)],[width3 heigth3]); 
gauss_filt_3 = imgaussfilt(rend3,sigma);
I32_3=uint32(gauss_filt_3);



figure('Position',[500 500 500 900])

subplot(3,2,1)
scatter(peaks1(:,xCol),peaks1(:,yCol),1,'black');
axis([min(peaks1(:,xCol)) max(peaks1(:,xCol)) min(peaks1(:,yCol)) max(peaks1(:,yCol))])
axis square
% axis([600 1.2e3 1.032e5 1.04e5])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,2)
imshow(imrotate(I32_1,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
%axis square
axis off

subplot(3,2,3)
scatter(peaks2(:,xCol),peaks2(:,yCol),1,'black');
axis([min(peaks2(:,xCol)) max(peaks2(:,xCol)) min(peaks2(:,yCol)) max(peaks2(:,yCol))])
axis square
% axis([600 1.2e3 1.032e5 1.04e5])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,4)
imshow(imrotate(I32_2,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

subplot(3,2,5)
scatter(peaks3(:,xCol),peaks3(:,yCol),1,'black');
axis([min(peaks3(:,xCol)) max(peaks3(:,xCol)) min(peaks3(:,yCol)) max(peaks3(:,yCol))])
axis square
% axis([600 1.2e3 1.032e5 1.04e5])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,6)
imshow(imrotate(I32_3,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off


%% Plot MT Data from measurement

cent1 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_MT.mat');
cent2 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_MT_2.mat');
cent3 = load('/Users/Christian/GitHub/SMLM_visualization/test_data/real_MT_large_FOV.mat');

% x,y,photons, uncertainty, frame

peaks1 = cent1.subset2;
peaks2 = cent2.subset2;
peaks3 = cent3.subset2;

xCol = 2;
yCol = 3;
frameCol = 5;
photonsCol = 7;

% Plot the raw data for Cent1

pxlsize=10;
sigma = 0.7; 

heigth1=round((max(peaks1(:,yCol))-min(peaks1(:,yCol)))/pxlsize);
width1=round((max(peaks1(:,xCol))-min(peaks1(:,xCol)))/pxlsize);

rend1 = hist3([peaks1(:,xCol),peaks1(:,yCol)],[width1 heigth1]); 
gauss_filt_1 = imgaussfilt(rend1,sigma);
I32_1=uint32(gauss_filt_1);

% Plot the raw data for Cent2

heigth2=round((max(peaks2(:,yCol))-min(peaks2(:,yCol)))/pxlsize);
width2=round((max(peaks2(:,xCol))-min(peaks2(:,xCol)))/pxlsize);

rend2 = hist3([peaks2(:,xCol),peaks2(:,yCol)],[width2 heigth2]); 
gauss_filt_2 = imgaussfilt(rend2,sigma);
I32_2=uint32(gauss_filt_2);

% Plot the raw data for Cent3

heigth3=round((max(peaks3(:,yCol))-min(peaks3(:,yCol)))/pxlsize);
width3=round((max(peaks3(:,xCol))-min(peaks3(:,xCol)))/pxlsize);

rend3 = hist3([peaks3(:,xCol),peaks3(:,yCol)],[width3 heigth3]); 
gauss_filt_3 = imgaussfilt(rend3,sigma);
I32_3=uint32(gauss_filt_3);



figure('Position',[500 500 500 900])

subplot(3,2,1)
scatter(peaks1(:,xCol)-min(peaks1(:,xCol)),peaks1(:,yCol)-min(peaks1(:,yCol)),1,'black');
% axis([min(peaks1(:,xCol))-min(peaks1(:,xCol)) max(peaks1(:,xCol))-min(peaks1(:,xCol)) min(peaks1(:,yCol))-min(peaks1(:,yCol)) max(peaks1(:,yCol))-min(peaks1(:,yCol))])
axis square
axis([0 2000 -100 400])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,2)
imshow(imrotate(I32_1,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
%axis square
axis off

subplot(3,2,3)
scatter(peaks2(:,xCol)-min(peaks2(:,xCol)),peaks2(:,yCol)-min(peaks2(:,yCol)),1,'black');
% axis([min(peaks2(:,xCol)) max(peaks2(:,xCol)) min(peaks2(:,yCol)) max(peaks2(:,yCol))])
axis square
axis([-200 1000 0 2000])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,4)
imshow(imrotate(I32_2,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

subplot(3,2,5)
scatter(peaks3(:,xCol)-min(peaks3(:,xCol)),peaks3(:,yCol)-min(peaks3(:,yCol)),1,'black');
% axis([min(peaks3(:,xCol)) max(peaks3(:,xCol)) min(peaks3(:,yCol)) max(peaks3(:,yCol))])
axis square
axis([0 3000 0 8000])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,6)
imshow(imrotate(I32_3,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off
