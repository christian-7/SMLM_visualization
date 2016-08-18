clc, clear, close all

%% Plot Data from SureSim

cd('.\simulated_test_data');

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

% Plot the raw data

pxlsize = 10;
sigma   = 0.7; 

heigth=round((max(peaks(:,yCol))-min(peaks(:,yCol)))/pxlsize);
width=round((max(peaks(:,xCol))-min(peaks(:,xCol)))/pxlsize);

rendered_image = hist3([peaks(:,xCol),peaks(:,yCol)],[width heigth]); 
gauss_filt = imgaussfilt(rendered_image,sigma);
I32=uint32(gauss_filt);

figure('Position',[500 500 800 500],'name','sureSim')
subplot(1,2,1)
scatter(peaks(:,xCol),peaks(:,yCol),1,'black.');
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

cd('..\')

%% Plot Data from SimpleSimulator

cd('.\simulated_test_data'); 

load('simulated_MT_3D_radius_20nm_015_wGT.mat')
peaks1 = sim_line;

load('simulated_MT_3D_radius_20nm_030_wGT.mat');
peaks2 = sim_line;

xCol = 1;
yCol = 2;
frameCol = 4;
photonsCol = 3;

% Plot the raw data

pxlsize = 5;
sigma   = 1.2 

heigth=round((max(peaks1(:,yCol))-min(peaks1(:,yCol)))/pxlsize);
width=round((max(peaks1(:,xCol))-min(peaks1(:,xCol)))/pxlsize);

rendered_image = hist3([peaks1(:,xCol),peaks1(:,yCol)],[width heigth]); 
gauss_filt = imgaussfilt(rendered_image,sigma);
I32=uint32(gauss_filt);

heigth=round((max(peaks2(:,yCol))-min(peaks2(:,yCol)))/pxlsize);
width=round((max(peaks2(:,xCol))-min(peaks2(:,xCol)))/pxlsize);

rendered_image2 = hist3([peaks2(:,xCol),peaks2(:,yCol)],[width heigth]); 
gauss_filt2 = imgaussfilt(rendered_image2,sigma);
I322=uint32(gauss_filt2);

figure('Position',[400 400 1200 300],'name','simpleSim')
subplot(2,2,1)
scatter(peaks1(:,xCol),peaks1(:,yCol),1,'black.');
axis([min(peaks1(:,xCol)) max(peaks1(:,xCol)) min(peaks1(:,yCol)) max(peaks1(:,yCol))])
% axis square
axis([0 1e3 -50 100])
box on
title('labelling efficiency = 0.15');
xlabel('nm');
ylabel('nm');

subplot(2,2,2)
imshow(imrotate(I32,90),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

subplot(2,2,3)
scatter(peaks2(:,xCol),peaks2(:,yCol),1,'black.');
axis([min(peaks2(:,xCol)) max(peaks2(:,xCol)) min(peaks2(:,yCol)) max(peaks2(:,yCol))])
% axis square
axis([0 1e3 -50 100])
box on
title('labelling efficiency = 0.30');
xlabel('nm');
ylabel('nm');

subplot(2,2,4)
imshow(imrotate(I322,90),[0 30]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

cd('..\');

%% Plot Centriole Data from measurement

cd('.\test_data'); 

cent1 = load('real_Cent_1.mat');
cent2 = load('real_Cent_2.mat');
cent3 = load('real_Cent_large.mat');

% x,y,photons, uncertainty, frame

peaks1 = cent1.subset2;
peaks2 = cent2.subset2;
peaks3 = cent3.subset2;

xCol = 1;
yCol = 2;
frameCol = 5;
photonsCol = 3;

% Plot the raw data for Cent1

pxlsize = 5;
sigma   = 1.2; 

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



figure('Position',[200 200 500 900],'name','real Centrioles')

subplot(3,2,1)
scatter(peaks1(:,xCol)-min(peaks1(:,xCol)),peaks1(:,yCol)-min(peaks1(:,yCol)),1,'black.');
axis square
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,2)
imshow(imrotate(I32_1,90),[0 10]);
colormap hot
title('Blurred 2D Histogram');
%axis square
axis off

subplot(3,2,3)
scatter(peaks2(:,xCol)-min(peaks2(:,xCol)),peaks2(:,yCol)-min(peaks2(:,yCol)),1,'black.');
axis square
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,4)
imshow(imrotate(I32_2,90),[0 10]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

subplot(3,2,5)
scatter(peaks3(:,xCol)-min(peaks3(:,xCol)),peaks3(:,yCol)-min(peaks3(:,yCol)),1,'black.');
axis square
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,6)
imshow(imrotate(I32_3,90),[0 10]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

cd('..\')

%% Plot MT Data from measurement

cd('.\test_data'); 

cent1 = load('real_MT.mat');
cent2 = load('real_MT_2.mat');
cent3 = load('real_MT_large_FOV.mat');

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



figure('Position',[500 100 1200 500],'name','real MT')

subplot(3,2,1)
scatter(peaks1(:,xCol)-min(peaks1(:,xCol)),peaks1(:,yCol)-min(peaks1(:,yCol)),1,'black.');
% axis([min(peaks1(:,xCol))-min(peaks1(:,xCol)) max(peaks1(:,xCol))-min(peaks1(:,xCol)) min(peaks1(:,yCol))-min(peaks1(:,yCol)) max(peaks1(:,yCol))-min(peaks1(:,yCol))])
% axis square
% axis([0 2000 -100 400])
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
scatter(peaks2(:,yCol)-min(peaks2(:,yCol)),peaks2(:,xCol)-min(peaks2(:,xCol)),1,'black.');
% axis([min(peaks2(:,xCol)) max(peaks2(:,xCol)) min(peaks2(:,yCol)) max(peaks2(:,yCol))])
% axis square
% axis([-200 1000 0 2000])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,4)
% imshow(imrotate(I32_2,360),[0 20]);
imshow(flipdim(I32_2,1),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

subplot(3,2,5)
scatter(peaks3(:,yCol)-min(peaks3(:,yCol)),peaks3(:,xCol)-min(peaks3(:,xCol)),1,'black.');
% axis([min(peaks3(:,xCol)) max(peaks3(:,xCol)) min(peaks3(:,yCol)) max(peaks3(:,yCol))])
% axis square
axis([0 8000 0 3000])
box on
title('Scatter plot');
xlabel('nm');
ylabel('nm');

subplot(3,2,6)
% imshow(imrotate(I32_3,90),[0 20]);
imshow(flipdim(I32_3,1),[0 20]);
colormap hot
title('Blurred 2D Histogram');
% axis square
axis off

cd('..\')