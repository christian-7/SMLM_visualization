clear, clc, close all

%% Generate randomly localized points along a line
% these are the true molecule positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length_of_line = 1000;
spacing = 5;
line1 = zeros(length_of_line/spacing,2);
line1(1,:)=[1,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(line1);

line1(i,1) = line1(i-1,1)+spacing;
    
end

line1(:,2) = 490;
line2 = line1;
line2(:,2) = 510;
mol_list = [line1;line2];

figure
scatter(mol_list(:,1),mol_list(:,2),1,'red');hold on;
axis([0 1000 0 1000])

%% Around each molecule generate localizations 

% Generate a uniform distribution between 0 and 1

uni_dis  = makedist('Uniform','lower',0,'upper',1);

% Compute the cdfs between 0 and 1, 
% cdf1 = list of 1e6 random numbers between 0and 1

x = 0:1e-6:1;
y = pdf(uni_dis,x);
cdf1 = cdf(uni_dis,x);

%% How many locs per molecule

rando_N_of_locs = [];

for i = 1:length(mol_list);
        
    r = [];
    r = randi([0,1e6],1,1);
    rando_N_of_locs(i,1) = cdf1(r);
    
end

%% Simulate locslizations around each molecule

% 1. On Time
load('Z:\Christian-Sieben\data_HTP\2016-04-13_A647_Cent_Calibration\locResults\NB_Nbr5_A647_COT_1500mW_1\photophysics_FOV1\locs_per_mol_1.mat')
dist_locs = fitdist(nbr_of_locs, 'lognormal');

% 2. Position x and y
load('Z:\Christian-Sieben\data_HTP\2016-04-13_A647_Cent_Calibration\locResults\NB_Nbr5_A647_COT_1500mW_1\photophysics_FOV1\radius1.mat')
dist_xpos = fitdist(allclustersCx, 'Normal');
dist_ypos = fitdist(allclustersCy, 'Normal');

% 2. Photons
load('Z:\Christian-Sieben\data_HTP\2016-04-13_A647_Cent_Calibration\locResults\NB_Nbr5_A647_COT_1500mW_1\photophysics_FOV1\photons1.mat')
dist_pho  = fitdist(pho,'Kernel','Width',100);
% x = 100:10:8000;
% y = pdf(dist_pho,x);
% scatter(x,y)

% 2. Dark Time
load('Z:\Christian-Sieben\data_HTP\2016-04-13_A647_Cent_Calibration\locResults\NB_Nbr5_A647_COT_1500mW_1\photophysics_FOV1\dt1_in_frames.mat')
dist_dT  = fitdist(allgaps,'Kernel','Width',100);
% x = 100:10:8000;
% y = pdf(dist_dT,x);
% scatter(x,y)

%% 

all_sim_x = [];
all_sim_y = [];
all_sim_photons = [];
all_sim_frame = [];

for i = 1:length(mol_list);           % for all molecules
    
    % x,y,photons,frame
    
    sim_locs = zeros(round(icdf(dist_locs,rando_N_of_locs(i,1))),4);            % generate an empty container for each molecule
    startFrame = randi([0,1e4],1,1);                                            % appearance of the first loc for the respective molecule
    
                        rando_N = [];                                           % generate random numbers to pick the parameters for each localization
     
                        for k = 1:length(sim_locs(:,1));
        
                        r = [];
                        r = randi([0,1e6],1,1);
                        rando_N(k,1) = cdf1(r);
    
                        r = [];
                        r = randi([0,1e6],1,1);
                        rando_N(k,2) = cdf1(r);
   
                        r = [];
                        r = randi([0,1e6],1,1);
                        rando_N(k,3) = cdf1(r);
    
                        r = [];
                        r = randi([0,1e6],1,1);
                        rando_N(k,4) = cdf1(r);

                        end
    
 for j = 1:length(sim_locs(:,1));                                       % for all locs     
    
    sim_locs(j,1) = mol_list(i,1) + icdf(dist_xpos,rando_N(j,1));       % x coordinate = true molecule position plus loc precision
    sim_locs(j,2) = mol_list(i,2) + icdf(dist_ypos,rando_N(j,2));       % y coordinate = true molecule position plus loc precision
    sim_locs(j,3) = icdf(dist_pho,rando_N(j,3));                        % photons
    sim_locs(j,4) = startFrame + round(icdf(dist_dT,rando_N(j,4)));     % start frame plus dT
    
 end         
    
            all_sim_x       = vertcat(all_sim_x,sim_locs(:,1));
            all_sim_y       = vertcat(all_sim_y,sim_locs(:,2));
            all_sim_photons = vertcat(all_sim_photons, sim_locs(:,3));
            all_sim_frame   = vertcat(all_sim_frame, sim_locs(:,4));
 
end

sim_line(:,1) = all_sim_x;
sim_line(:,2) = all_sim_y;
sim_line(:,3) = all_sim_photons;
sim_line(:,4) = all_sim_frame;


% x2 = 0:0.01:1;
% icdf1 = icdf(dist_xpos,x2);
% figure; scatter(x2,icdf1)

% pd2  = fitdist(tracklength,'Kernel','Width',0.2);
% x = 0:1:100;
% y = pdf(dist_pho,x);


%% 
 
scatter(mol_list(:,1),mol_list(:,2),1,'red');hold on;
scatter(all_sim_x,all_sim_y,1,'blue')
axis([0 1e3 200 600]);


%% Render the image

pxlsize=5;

heigth = round((max(all_sim_y)-min(all_sim_y))/pxlsize);
width  = round((max(all_sim_x)-min(all_sim_x))/pxlsize);
       
rendered = hist3([all_sim_y,all_sim_x],[height width]);

imshow(rendered);
colormap('hot');

%% 


name = ['image_5nm.tiff'];

I32=[];
I32=uint32(rendered);

t = Tiff(name,'w');

tagstruct.ImageLength       = size(I32,1);
tagstruct.ImageWidth        = size(I32,2);
tagstruct.Photometric       = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample     = 32;
tagstruct.SamplesPerPixel   = 1;
tagstruct.RowsPerStrip      = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()
