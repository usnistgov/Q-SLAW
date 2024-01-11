function mean_resid = compute_localization_residual(thresh)
% Compute the localization residual in a least squares sense from the
% actual onset to the detected onset.

%init vars
mc_runs = 30;
loc_sim_params = repmat(struct('sizeE',[1,1],'num_steps',1,'nu',1,'strain_step',1,...
    'thresh',1,'num_loc_bands',1,'localization_width',1,'onset',1,'L',1,'k',1),mc_runs,1);
loc_det_results = repmat(struct('onset_actual',1,'localization_str',1,'onset_pca',1),mc_runs,1);
onset_resid = zeros(1,mc_runs);

%run the MC for loop
for mc_count = 1:mc_runs
    
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
    onset_resid(mc_count) = (loc_det_results(mc_count).onset_wav - loc_det_results(mc_count).onset_actual).^2;

    if mod(mc_count,100) == 0
        disp(mc_count)
    end
    
end

mean_resid = mean(onset_resid);
