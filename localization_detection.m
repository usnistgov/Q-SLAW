function [loc_params] = localization_detection(loc_params,E_loc)
%
% Localization routine with wavelet transforms.
%
% For a MC approach simulation parameters are taken as input, default is:
%    %experiment setup params
%    loc_sim_params.sizeE = [128,512]; %size of strain field to simulate
%    loc_sim_params.num_steps = 150; %nubmer of total images
%    loc_sim_params.nu = 0.15; %effective linear elastic-like poisson's ratio
%    loc_sim_params.strain_step = 0.005; %strain per image
%
%    %localization band params
%    loc_sim_params.num_loc_bands = 12; % number of localization areas to use
%    loc_sim_params.localization_width = 2.0; % width of localization band
%                % (strength of sigma for imgaussfilt: = 1 => ~7px width)
%    loc_sim_params.onset = 0.7; %normalized location of onset
%    loc_sim_params.L = 0.3; %strength of localization bands
%    loc_sim_params.k = 0.15; %steepness of localization onset
%    loc_sim_params.thresh = 0.01; %strength threshold for localization onset
%
%   E_loc: data to run the algorithm on, leave blank to computationally 
%          generate a map with the loc_sim_params
%
%  OUTPUTS:
%  loc_params struct containing:
%   .localization_str: localization strength metric per strain step
%   .localization_str_uncert: localization strength metric uncerainty per strain step
%   .localization_amnt: number of localizations events per strain step
%   .localization_onset: estimated localation where the localization begins,
%                      usually when amnt drops and str increases
%
% Alex Landauer, NIST MML, MMSD, Nov 2023
%

%% Set defaults
if nargin < 1
    %experiment setup params
    loc_params.sizeE = [128,512]; %size of strain field to simulate
    loc_params.num_steps = 150; %number of total images
    loc_params.nu = 0.15; %effective linear elastic-like poisson's ratio
    loc_params.strain_step = 0.005; %strain per image

    %localization band params
    loc_params.num_loc_bands = 12; % number of localization areas to use
    loc_params.localization_width = 2.0; % width of localization band
    % (strength of sigma for imgaussfilt: = 1 => ~7px width)
    loc_params.onset = 0.7;  %normalized location of onset
    loc_params.L = 0.3; %strength of localization bands
    loc_params.k = 0.15; %steepness of localization onset
    
    loc_params.thresh = 0.75; %strength threshold for localization onset

end

%% generate a strain map if no data is given
if nargin < 2
    % Generate strain maps with known bands and band strength
    E_loc = gen_synth_localized_strain_maps(loc_params);
end


%%
% figure
% hold on
% for ii = 1:length(E_loc_sim)
%     imagesc(E_loc_sim{ii}{1,1}),axis image,colorbar
%     drawnow
% end

%% Compute the wavelet transform with a isotropic wavelet, and compute localization metrics

for  ii = 1:length(E_loc)
    %get nan mask area 
    nan_area = double(~isnan(E_loc{ii}{2,2}));
    nan_area(nan_area == 0) = nan;

    %remove NaNs and replace with mean value
    E_loc{ii}{2,2}(isnan(E_loc{ii}{2,2})) = mean(E_loc{ii}{2,2},'all','omitnan');
    
    %run the wavelet analysis
    wscales = [2,4,8]/(round(128/size(E_loc{1}{2,2},1)));
    if wscales(1) < 1
        wscales = [1,2,4];
    end
    cwtstruct{ii} = cwtft2(E_loc{ii}{2,2},wavelet="marr",norm="L2",scales=wscales);%,angles=[0 pi/2 pi]
    
    %two measures of localization - strength (dynamic range of the wavelet
    %response) and amount (number of distinct regions in the wavelet response) 
    se = strel('disk',1);
    for scls = 1:length(wscales)
        %for each scale, binarize 
        BW_wav = imerode(imbinarize((-cwtstruct{ii}.cfs(:,:,1,scls))/max(abs(cwtstruct{ii}.cfs(:,:,1,scls)),[],'all','omitnan')),se);

        %find the number of "large" regions
        BWCC = bwconncomp(nan_area.*BW_wav);
        lcl_amnt_(scls) = sum(cellfun(@length,BWCC.PixelIdxList) > 0.005*prod(BWCC.ImageSize));
        
        %find the difference in the mean strain inside the localization
        %regions and outside of the localization
        try
            if wscales(scls) > 2
                px_list = BWCC.PixelIdxList{cellfun(@length,BWCC.PixelIdxList) <= 0.005*prod(BWCC.ImageSize)};
            end
            BW_wav_d = double(BW_wav);
            BW_wav_d(px_list) = 0;
        catch
            BW_wav_d = double(BW_wav);
        end
        BW_wav_d(BW_wav_d == 0) = nan;
        BW_wav_inv_d = double(~imdilate(BW_wav,se));
        BW_wav_inv_d(BW_wav_inv_d == 0) = nan;

        %define the localization strength as the mean strain in the areas
        %id'd in the wavelet analysis minus the main strain + noise (1 std
        %dev) everywhere else.
        localization_str_(scls) = abs( mean(abs(nan_area.*BW_wav_d.*E_loc{ii}{2,2}),'all','omitnan') -...
            (abs(mean(nan_area.*BW_wav_inv_d.*E_loc{ii}{2,2},'all','omitnan'))));
        localization_str_uncert_(scls) = sqrt(std(abs(nan_area.*BW_wav_d.*E_loc{ii}{2,2}),[],'all','omitnan')^2 + ...
            std(nan_area.*BW_wav_inv_d.*E_loc{ii}{2,2},[],'all','omitnan')^2);

    end
    
    % use the localization amount in the first scale since defuse
    % localization occurs at higher freqs
    loc_params.localization_amnt(ii) = lcl_amnt_(1);%mean(lcl_amnt_);
    %Save the strength and uncertainty (from variablity of the fields)
    [loc_params.localization_str(ii),idx] = max(localization_str_,[],'omitnan');
    loc_params.localization_str_uncert(ii) = localization_str_uncert_(idx);
        
end

%%
% 
% figure
% sca=1;
% for ii = 1:length(E_loc)
%     % imagesc(abs(sum(cwtstruct.cfs(:,:,1,1:4,sca),4)))
%     % imagesc(cwtstruct{ii}.cfs),colorbar
%     imagesc(latent_params(:,:,ii)),colorbar
%     drawnow
% end

%% Find the "localization onset" 

%one metric - where the number of regions decreases
num_loc_drop_idx = findchangepts(loc_params.localization_amnt,MaxNumChanges=1);

%another metric - when strength is over a threshold - define this as 3x the
%background noise at biggest wavelet scale in strain at max strain or a predefined 
%threshold, whichever is bigger
loc_params.thresh = max(3*std(nan_area.*BW_wav_inv_d.*E_loc{ii}{2,2},[],'all','omitnan'),loc_params.thresh);
[idx_thresh] = find(loc_params.localization_str > loc_params.thresh,1);
% idx_thresh = findchangepts(localization_str,MaxNumChanges=1); % alternate method, doesn't work as well

% save the index at which localization is taken to have started from
if isempty(idx_thresh) 
    %no localization if there's no drop or there's no strength above the
    %threshold
    loc_params.localization_onset = loc_params.num_steps+200;
    disp('below thresh')

elseif isempty(num_loc_drop_idx)
    %no localization if there's no drop or there's no strength above the
    %threshold
    loc_params.localization_onset = loc_params.num_steps+200;
    disp('no drop')

else
    %if both metric are available, estimate the onset as the mean of the two
    loc_params.localization_onset = (idx_thresh+num_loc_drop_idx)/2;%mean([idx_thresh,num_loc_drop_idx]);

end

% figure
% plot(localization_str)

% figure
% imagesc(E_loc_sim{ii}{1,1})
end