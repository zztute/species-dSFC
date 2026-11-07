% RUN_COMPUTE_ROI_COORDINATES Select a representative voxel for each Schaefer ROI.
% Written by Junrui Zhang, 2025
%
% Inputs: Subject-specific Schaefer 400 atlases in diffusion-image space.
% Outputs: min_points.mat with min_points [400 x 3] and label_set [1 x 400].
% Method: Select the voxel minimizing summed within-ROI Euclidean distance.
% Coordinates: Voxel indices, following the supplied algorithm.
% Dependencies: SPM (spm_vol and spm_read_vols).
% Usage: Set data_root below, then run this script.

clear;
clc;
tic;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
sc_root = fullfile(data_root, 'SC', 'Subj'); % Subject folders containing structural connectomes and diffusion-space atlases.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the dSFC module.
n_nodes = 400;  % the number of parcels
label_set = 1:n_nodes;

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
assert(exist('spm_vol', 'file') == 2 && exist('spm_read_vols', 'file') == 2, ...
    'Add SPM to the MATLAB path; see README.md.');
subject_ids = get_subject_ids(sc_root);

%% Calculate ROI representative coordinates
for subject = 1:numel(subject_ids)
    subject_id = subject_ids{subject};
    fprintf('%s is processing...\n', subject_id);
    subject_dir = fullfile(sc_root, subject_id); % Structural-input and output directory for this subject.
    atlas_file = fullfile(subject_dir, 'atlas_in_b0_int.nii'); % Schaefer 400 atlas resampled to this subject's diffusion-image space.
    compressed_atlas = fullfile(subject_dir, 'atlas_in_b0_int.nii.gz'); % Optional compressed version of the diffusion-space atlas.
    if ~isfile(atlas_file)
        assert(isfile(compressed_atlas), 'Missing diffusion-space atlas for subject %s.', subject_id);
        gunzip(compressed_atlas, subject_dir);
    end
    volume = spm_vol(atlas_file);
    assert(isscalar(volume), 'The atlas must be a single three-dimensional volume.');
    atlas_data = spm_read_vols(volume);
    assert(all(isfinite(atlas_data(:))) && ...
        all(abs(atlas_data(:) - round(atlas_data(:))) < 1e-6), ...
        'Use an integer-label atlas resampled with nearest-neighbor interpolation.');
    roi_labels = round(atlas_data(:));
    assert(all(ismember(unique(roi_labels), [0, label_set])), ...
        'The diffusion-space atlas may contain only labels 0-400.');
    [x, y, z] = ndgrid(1:volume.dim(1), 1:volume.dim(2), 1:volume.dim(3));
    voxel_coordinates = [x(:), y(:), z(:)];
    min_points = zeros(n_nodes, 3);
    for region = 1:n_nodes
        selected_coordinates = voxel_coordinates(roi_labels == label_set(region), :);
        assert(~isempty(selected_coordinates), ...
            'ROI %d has no voxels for subject %s.', label_set(region), subject_id);
        min_points(region, :) = compute_roi_medoid(selected_coordinates);
    end
    output_file = fullfile(subject_dir, 'min_points.mat'); % Representative ROI voxel coordinates and their label order.
    coordinate_units = 'voxel_indices';
    save(output_file, 'min_points', 'label_set', 'coordinate_units');
    fprintf('%s is done.\n', subject_id);
end
toc;
