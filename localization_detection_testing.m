
%
% Localization detection testing with Monte Carlo
%
% randomize the simulation parameters and check that onset and strength
% remain accurate
%
%
% Alex Landauer, NIST MML, MMSD, Sept 2023
%

clear, close all

%Monte Carlo setup params
rng(4);
mc_runs = 300;

%% find the best threshold

%use fminsearch to find the best threshold value, and plot the objective
%func values
options = optimset('PlotFcns',@optimplotfval);
thresh = fminsearch(@compute_localization_residual,0.048,options);
% thresh = 0.0486;

%% Run the MC for loop with that threshold
loc_sim_params = repmat(struct('sizeE',[1,1],'num_steps',1,'nu',1,'strain_step',1,...
    'thresh',1,'num_loc_bands',1,'localization_width',1,'onset',1,'L',1,'k',1),mc_runs,1);
loc_det_results = repmat(struct('onset_actual',1,'localization_str',1,'onset_pca',1),mc_runs,1);
onset_resid = zeros(1,mc_runs);

parfor mc_count = 1:mc_runs
    
    %experiment setup params
    loc_sim_params(mc_count).sizeE = [128,512]; %size of strain field to simulate
    loc_sim_params(mc_count).num_steps = 130; %number of total images
    loc_sim_params(mc_count).nu = 0.15; %effective linear elastic-like poisson's ratio
    loc_sim_params(mc_count).strain_step = 0.005; %strain per image

    loc_sim_params(mc_count).thresh = thresh; %localization onset strength threshold 

    %localization band params
    loc_sim_params(mc_count).num_loc_bands = randi([1,12]); % number of localization areas to use
    loc_sim_params(mc_count).localization_width = randi([1,50])/4; % width of localization band
    loc_sim_params(mc_count).onset = randi([15,80])/100;  %normalized location of onset
    loc_sim_params(mc_count).L = randi([16,60])/80; %strength of localization bands
    loc_sim_params(mc_count).k = randi([5,35])/100; %steepness of localization onset

    %run the localization detection
    loc_det_results(mc_count).onset_actual = loc_sim_params(mc_count).num_steps*loc_sim_params(mc_count).onset;
    loc_params = localization_detection(loc_sim_params(mc_count));
    loc_det_results(mc_count).localization_str = loc_params.localization_str;
    loc_det_results(mc_count).localization_amnt = loc_params.localization_amnt;
    loc_det_results(mc_count).onset_wav = loc_params.localization_onset;

    %compute the residual from the detected value and the nominal
    onset_resid(mc_count) = (loc_det_results(mc_count).onset_wav - loc_det_results(mc_count).onset_actual);

    % if mod(mc_count,100) == 0
    %     disp(mc_count)
    % end
    
end

%% plot results
figure
histogram(onset_resid,20)
% mean(onset_resid)
% std(onset_resid)

for ii = 1:mc_runs
    loc_str_map(:,ii) = loc_det_results(ii).localization_str;
    loc_amt_map(:,ii) = loc_det_results(ii).localization_amnt;

    STR_mea(ii) = loc_det_results(ii).localization_str(end);
    STR_sim(ii) = loc_sim_params(ii).L;
    AMNT_sim(ii) = loc_sim_params(ii).num_loc_bands;
end

figure
imagesc(loc_str_map'),colorbar

figure
imagesc(loc_amt_map'),colorbar

figure
plot(AMNT_sim)
hold on
plot(loc_amt_map(end,:))

onset_err_mean = mean(onset_resid(onset_resid<100));
onset_err_std = std(onset_resid(onset_resid<100));
