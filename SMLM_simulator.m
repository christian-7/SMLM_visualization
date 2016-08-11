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

%% Generate a 3D projection of a microtubule
% these are the true molecule positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angle = 1:(360/13):360; % 13 monomers per circle, divide the circle in 13 parts (angles)
radius = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = sind(90-angle)*radius; 
c = radius-b;
c = c(:);

% % Plot the distribution for one circle
% c(:,2)=1;
% scatter(c(:,1),c(:,2));

% Multiply rings n-times
m = 1;
n = 4;
o = 13;

mol_list = [];
mol_list(1:13,1) = 1;
mol_list(1:13,2) = c(1:13,1);

for i = 1:250; % i*4 = length of the simulated MT i nm

mol_list((i*o)+1:(i*o)+o,1) = i*n;
mol_list((i*o)+1:(i*o)+o,2) = mol_list(1:13,2);

end

figure
scatter(mol_list(:,1),mol_list(:,2),5,'*r');
axis([0 1000 -50 100]);
box on;
axis square;

%% Reduce the amount of true molecules to account for labeling efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelling_eff = 0.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbr_of_labels = length(mol_list)*labelling_eff;

% Select random molecules fro the list

r = randi([1,length(mol_list)],1,round(nbr_of_labels));

mol_list2 = [];
mol_list2 = mol_list(r,1:end);

figure
scatter(mol_list(:,1),mol_list(:,2),5,'*r'); hold on;
scatter(mol_list2(:,1),mol_list2(:,2),5,'ob');
axis([0 1000 -50 100]);
legend('all','subset');
box on;
axis square

%% Generate a uniform distribution between 0 and 1

uni_dis  = makedist('Uniform','lower',0,'upper',1);

% Compute the cdfs between 0 and 1, 
% cdf1 = list of 1e6 random numbers between 0and 1

x = 0:1e-6:1; 
y = pdf(uni_dis,x);
cdf1 = cdf(uni_dis,x); % generate 1e6 random numbers between 0 and 1, this is later used to find random points on the experimental icdf

%% Simulate localizations around each molecule

% 1. On Time
load('K:\Christian\GitHub\SMLM_vis\exp_dist\locs_per_mol_1.mat')
dist_locs = fitdist(nbr_of_locs, 'lognormal');

% 2. Position x and y
load('K:\Christian\GitHub\SMLM_vis\exp_dist\radius1.mat')
dist_xpos = fitdist(allclustersCx, 'Normal');
dist_ypos = fitdist(allclustersCy, 'Normal');

% 2. Photons
load('K:\Christian\GitHub\SMLM_vis\exp_dist\photons1.mat')
dist_pho  = fitdist(pho,'Kernel','Width',100);
% x = 100:10:8000;
% y = pdf(dist_pho,x);
% scatter(x,y)

% 2. Dark Time
load('K:\Christian\GitHub\SMLM_vis\exp_dist\dt1_in_frames.mat')
dist_dT  = fitdist(allgaps,'Kernel','Width',100);
% x = 100:10:8000;
% y = pdf(dist_dT,x);
% scatter(x,y)

%% Start Simulator

nframes = 1e4;

% 1. How many locs per molecule

% for each true molecule, it picks a random number 
% between 0 and 1 from the uniform CDF

rando_N_of_locs = [];

for i = 1:length(mol_list2);
        
    r = [];
    r = randi([0,1e6],1,1);
    rando_N_of_locs(i,1) = cdf1(r);
    
end

tic

rng('shuffle')

all_sim_x = [];
all_sim_y = [];
all_sim_photons = [];
all_sim_frame = [];

% 2. Simulate localizations for each true molecule 

for i = 1:length(mol_list2);           % for all molecules
    
    % x,y,photons,frame
    
    sim_locs = zeros(round(icdf(dist_locs,rando_N_of_locs(i,1))),4);            % generate an empty container for each molecule
    startFrame = randi([1,nframes],1,1);                                        % appearance of the first loc for the respective molecule, between 1- frames
    
                        rando_N = [];                                           % generate random numbers to pick the parameters for each localization
     
                        for k = 1:length(sim_locs(:,1));
        
                        r = [];
                        r = randi([0,1e6],1,4);
                        rando_N(k,1:4) = cdf1(r);

                        end
    
 for j = 1:length(sim_locs(:,1));                                        % for all locs     
    
    sim_locs(j,1) = mol_list2(i,1) + icdf(dist_xpos,rando_N(j,1));       % x coordinate = true molecule position plus loc precision
    sim_locs(j,2) = mol_list2(i,2) + icdf(dist_ypos,rando_N(j,2));       % y coordinate = true molecule position plus loc precision
    sim_locs(j,3) = icdf(dist_pho,rando_N(j,3));                         % photons
    sim_locs(j,4) = startFrame + round(icdf(dist_dT,rando_N(j,4)));      % start frame plus dT
    
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

fprintf(' -- Simulation done after %f sec -- \n',toc)

%% Plot the simulation result and the ground truth

figure 
scatter(mol_list2(:,1),mol_list2(:,2),1,'red');hold on;
scatter(all_sim_x,all_sim_y,1,'blue')
axis([0 1000 -50 50]);


%% Render the image

pxlsize=10;

heigth = round((max(all_sim_y)-min(all_sim_y))/pxlsize);
width  = round((max(all_sim_x)-min(all_sim_x))/pxlsize);
       
rendered = hist3([all_sim_y,all_sim_x],[heigth width]);

imshow(rendered);
colormap('hot');

%% Save image as 32-bit image

name = ['image_' num2str(pxlsize) 'nm_p_pxl.tiff'];

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

