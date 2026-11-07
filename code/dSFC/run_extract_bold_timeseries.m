% RUN_EXTRACT_BOLD_TIMESERIES Extract regional mean BOLD time series using Schaefer 400.
% Written by Junrui Zhang, 2025
%
% Inputs: Preprocessed four-dimensional fMRI images and a matching atlas.
% Outputs: One MAT file per subject, containing timeseries [time x 400].
% Dependencies: SPM (spm_vol and spm_read_vols).
% Usage: Set data_root below, then run this script.

clear;
clc;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
fmri_root = fullfile(data_root, 'BOLD', 'Subj'); % Subject folders containing preprocessed fMRI images.
atlas_file = fullfile(data_root, 'atlas', 'Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii'); % Schaefer 400 atlas in the fMRI image space.
output_dir = fullfile(data_root, 'BOLD', 'timeseries'); % Output directory for regional BOLD time-series MAT files.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the dSFC module.
n_nodes = 400;  % the number of parcels
label_set = 1:n_nodes;

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
assert(exist('spm_vol', 'file') == 2 && exist('spm_read_vols', 'file') == 2, ...
    'Add SPM to the MATLAB path; see README.md.');
assert(isfile(atlas_file), 'The Schaefer 400 atlas file is missing.');
subject_ids = get_subject_ids(fmri_root);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

%% Prepare atlas masks
V_atlas = spm_vol(atlas_file);
assert(isscalar(V_atlas), 'The atlas must be a single three-dimensional volume.');
atlas_data = spm_read_vols(V_atlas);
assert(all(isfinite(atlas_data(:))) && ...
    all(abs(atlas_data(:) - round(atlas_data(:))) < 1e-6), ...
    'The atlas must contain finite integer labels; use nearest-neighbor resampling.');
assert(isequal(setdiff(unique(atlas_data(:)), 0), label_set'), ...
    'The atlas must contain labels 1-400 and background label 0.');
masks = cell(n_nodes, 1);
for region = 1:n_nodes
    masks{region} = find(atlas_data == label_set(region));
end

%% Extract time series
for subject = 1:numel(subject_ids)
    subject_id = subject_ids{subject};
    fprintf('%s is processing...\n', subject_id);
    fmri_file = fullfile(fmri_root, subject_id, 'rest_proc_network.nii'); % Preprocessed fMRI image for this subject.
    assert(isfile(fmri_file), 'Missing fMRI input for subject %s.', subject_id);
    V_fmri = spm_vol(fmri_file);
    for timepoint = 1:numel(V_fmri)
        assert(isequal(V_atlas.dim, V_fmri(timepoint).dim) && ...
            max(abs(V_atlas.mat(:) - reshape(V_fmri(timepoint).mat, [], 1))) < 1e-4, ...
            'Atlas and fMRI images must have matching dimensions and spatial transforms.');
    end
    fmri_data = spm_read_vols(V_fmri);
    n_timepoints = numel(V_fmri);
    timeseries = zeros(n_timepoints, n_nodes);
    for timepoint = 1:n_timepoints
        volume = fmri_data(:, :, :, timepoint);
        for region = 1:n_nodes
            timeseries(timepoint, region) = mean(volume(masks{region}));
        end
    end
    assert(all(isfinite(timeseries(:))), 'Nonfinite BOLD values found for subject %s.', subject_id);
    output_file = fullfile(output_dir, [subject_id, '_BOLD_timeseries.mat']); % Regional mean BOLD time series, with ROI labels in label_set.
    save(output_file, 'timeseries', 'label_set');
    fprintf('%s is done.\n', subject_id);
end

% Next step: calculate dynamic FC from these time series using DynamicBC.
% This repository reads the resulting FCM structure; see README.md.
