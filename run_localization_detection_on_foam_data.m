%
% Localization analysis for example foam data
%
%
% Alex Landauer, NIST MML, MMSD, Jan 2024
%

%% set up inputs and algorithm settings
clear, close all

% ***** User settings: data location, intensity threshold, and filter type ***** 

% set up inputs
data_dir = 'example_dic_data';
qDIC_data = dir(['.',filesep,data_dir,filesep,'*run_data.mat']);

% algorithm settings
thresh = 0.05;
strain_filter = 'optimal9';

% ***** 

for ii = 1:length(qDIC_data)
    %% read in data
    %complete current dataset
    cur_data = load([data_dir,filesep,qDIC_data(ii).name]);
    cur_data.spc_prefix = qDIC_data(ii).name(1:end-13);

    %truncate to the first half (localization that we care about only occurs in
    %the loading segment
    u_curr = cur_data.u(1:ceil(length(cur_data.u)/2));
    %can import time if needed (e.g. for plotting)
    %time = cur_data.img_times(1:ceil(length(cur_data.u)/2));
    
    %Basic noise floor from first frame of displacement
    noise_floor = std(u_curr{1}{2},[],"all","omitmissing");

    %% Compute strain fields

    %cumulate displacement (maybe not needed?)
    [u_total,~] = inc2cum(u_curr,cur_data.dm,cur_data.gridPoints,'cubic');

    % compute strains
    Fij_total = calculateFij_2d(u_total,cur_data.dm,cur_data.scale,strain_filter);
    [Eij_total,eij_total] = calculateEij_2d(Fij_total);

    %compute medians
    for step = 1:length(Eij_total)
        med_strain22_all{ii}(step) = median(eij_total{step}{2,2},'all','omitmissing');
        std_strain22_all{ii}(step) = std(eij_total{step}{2,2},[],'all','omitmissing');
    end


    %% run the localization detection
    %experiment setup params, comment for an MC run
    loc_params{ii}.sizeE = size(Eij_total{1}{2,2}); %size of strain field to simulate
    loc_params{ii}.num_steps = length(Eij_total); %number of total images
    
    % MC simulation params, uncomment for MC run
    %     loc_sim_params(mc_count).nu = 0.15; %effective linear elastic-like poisson's ratio
    %     loc_sim_params(mc_count).strain_step = 0.005; %strain per image

    loc_params{ii}.thresh = thresh; %localization onset total strength threshold

    % MC sim localization band params, uncomment for MC run
    % loc_sim_params(mc_count).num_loc_bands = randi([1,40]); % number of localization areas to use
    % loc_sim_params(mc_count).localization_width = randi([1,50])/4; % width of localization band
    % loc_sim_params(mc_count).onset = randi([10,85])/100;  %normalized location of onset
    % loc_sim_params(mc_count).L = randi([1,25])/40; %strength of localization bands
    % loc_sim_params(mc_count).k = randi([5,35])/100; %steepness of localization onset

    loc_params{ii} = localization_detection(loc_params{ii},Eij_total);


    %plot the results
    figure(1)
    hold on

    [~,max_idx] = max(abs(med_strain22_all{ii})); %find the max strain and truncate there
    errorbar(-med_strain22_all{ii}(1:max_idx),loc_params{ii}.localization_str(1:max_idx),...
        loc_params{ii}.localization_str_uncert(1:max_idx),'Marker','*')
    xlabel('Median axial engineering strain, e_1_1')
    ylabel('localization strength metric')
    legend_labels1{ii} = cur_data.spc_prefix;
    

    figure(2)
    hold on
    plot(-med_strain22_all{ii}(1:max_idx),loc_params{ii}.localization_amnt(1:max_idx))
    xlabel('step')
    ylabel('localization amount metric')
    legend_labels2{ii} = cur_data.spc_prefix;
    legend(legend_labels2{ii},"Interpreter","none")

end
figure(1),legend(legend_labels1,"Interpreter","none")
figure(2),legend(legend_labels2,"Interpreter","none")
