function E_loc_sim = gen_synth_localized_strain_maps(params)
%generate synthetic data with strain localization strength increasing as
% a logistic function of number of steps
% 
% INPUT (as a parameter structure):
% sizeE: size of the strain map
% num_steps: number of strain steps to simulate (0.005 global strain per step)
% num_loc_bands: number of localization "bands" in the map
% localization_width: width of localization band (strength of sigma for
%                     imgaussfilt, so = 1 => ~7px width)
% nu: effective linear elastic-like poisson's ratio
% strain_step: strain per image
% L: strength of localization bands
% k: steepness of localization onset
%
%
% OUTPUT:
% E_loc_sim: simulated strain field with localization
%
% -- Alex Landauer, NIST MML, 2023-06-08
%


% ------ parameters ------
sizeE = params.sizeE;
num_steps = params.num_steps;
num_loc_bands = params.num_loc_bands;
localization_width = params.localization_width;
strain_step  = params.strain_step;
nu = params.nu;
L = params.L;
k = params.k;
onset = params.onset;
x0 = round(onset*num_steps); %strain step for midpoint of localization


% ------ localization domain ------

%make bands of semi-randomized width and location throughout the strain maps,
% which are constrained to limit overlap
band_domains = zeros(sizeE);
band_loc_centers = randperm(sizeE(1),num_loc_bands);
band_x_start = randi(sizeE(2)-round(sizeE(2)/8),[num_loc_bands,1]);
for ii = 1:num_loc_bands
    band_width = randi(round(sizeE(2) - band_x_start(ii)-1),[num_loc_bands,1]);
    band_domains(band_loc_centers(ii),band_x_start(ii):(band_x_start(ii)+band_width(ii))) = randi([1,20])/20;
end

%Filter the bands - this makes the localization regime be
%wider than a single px with a Gaussian decay. Normalize afterword to keep
%it in the magnitudes we want
band_domains = normalize(imgaussfilt(band_domains,localization_width),'range');


% ------ create the localized strain maps ------
ideal_strain_tensor = cell(1,num_steps);
localization_strain_tensor = cell(1,num_steps);
E_loc_sim = cell(1,num_steps);
for time_step = 1:num_steps
    %generate the ideal (uniform) strain map for the current strain step
    %with some white noise
    ideal_strain_tensor{time_step}{2,2} = imnoise((strain_step*time_step)*ones(sizeE),'speckle',0.001);
    ideal_strain_tensor{time_step}{1,2} = imnoise(zeros(sizeE),'speckle',0.001);
    ideal_strain_tensor{time_step}{2,1} = imnoise(zeros(sizeE),'speckle',0.001);
    ideal_strain_tensor{time_step}{1,1} = imnoise((nu*strain_step*time_step)*ones(sizeE),'speckle',0.001);
    

    %generate the localization band at the determined locations, which
    %increase in magnitude (stength) according to a logistic function of
    %the strain step
    localization_strain_tensor{time_step}{2,2} = ...
        (strain_step*time_step*localization_str(time_step,L,k,x0)).*band_domains;
    localization_strain_tensor{time_step}{1,2} = zeros(sizeE);
    localization_strain_tensor{time_step}{2,1} = zeros(sizeE);
    localization_strain_tensor{time_step}{1,1} = zeros(sizeE);
    
    %add the localizations to the ideal (uniform) strain map
    for ii = 1:2
        for jj = 1:2
            E_loc_sim{time_step}{ii,jj} = ideal_strain_tensor{time_step}{ii,jj} + ...
                localization_strain_tensor{time_step}{ii,jj};
        end
    end
end

end

function str = localization_str(time,L,k,x0)
%logistic function parameterized by inputs
% L = logistic strength
% k = sharpness
% x0 = centerpoint
%see: https://en.wikipedia.org/wiki/Logistic_function
str = L./(1+exp(-k*(time - x0)));

end


