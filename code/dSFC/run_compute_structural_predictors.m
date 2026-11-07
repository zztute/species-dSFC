% RUN_COMPUTE_STRUCTURAL_PREDICTORS Calculate structural-network predictors for Schaefer 400.
% Written by Junrui Zhang, 2025
%
% Inputs: FN_connectome_sift2.csv [400 x 400] and min_points.mat [400 x 3].
% Outputs: SC.mat, SC_bi.mat, pathlength.mat, communicability.mat,
%          and euclidean_dis.mat in each structural subject folder.
% Method: Binary shortest paths, binary communicability, and ROI distances.
% Dependencies: Brain Connectivity Toolbox; Statistics and Machine Learning Toolbox.
% Usage: Set data_root below; run run_compute_roi_coordinates.m first.

clear;
clc;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
sc_root = fullfile(data_root, 'SC', 'Subj'); % Subject folders containing structural-connectome CSV files and ROI coordinates.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the dSFC module.
n_nodes = 400;  % the number of parcels
label_set = 1:n_nodes;

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
assert(exist('distance_bin', 'file') == 2, 'Add Brain Connectivity Toolbox to the MATLAB path.');
subject_ids = get_subject_ids(sc_root);

%% Calculate structural predictors
for subject = 1:numel(subject_ids)
    subject_id = subject_ids{subject};
    fprintf('%s is processing...\n', subject_id);
    subject_dir = fullfile(sc_root, subject_id); % Structural-input and output directory for this subject.
    connectome_file = fullfile(subject_dir, 'FN_connectome_sift2.csv'); % Weighted structural connectome in ROI label order 1-400.
    coordinates_file = fullfile(subject_dir, 'min_points.mat'); % Representative ROI coordinates produced by run_compute_roi_coordinates.m.
    assert(isfile(connectome_file) && isfile(coordinates_file), ...
        'Missing structural connectome or ROI coordinates for subject %s.', subject_id);
    raw_connectome = readmatrix(connectome_file);
    assert(isequal(size(raw_connectome), [n_nodes, n_nodes]), ...
        'The structural connectome must be a 400-by-400 matrix.');
    assert(all(isfinite(raw_connectome(:))) && all(raw_connectome(:) >= 0), ...
        'Structural-connectome weights must be finite and nonnegative.');
    assert(norm(raw_connectome - raw_connectome', 'fro') <= ...
        1e-8 * max(1, norm(raw_connectome, 'fro')), ...
        'The input structural connectome must be symmetric.');
    SC_L2R_order = raw_connectome(label_set, label_set);
    SC_bi = double(SC_L2R_order > 0);
    sp = pathlength(SC_bi);
    co = fcn_communicability(SC_bi);

    coordinates = load(coordinates_file, 'min_points', 'label_set');
    assert(isfield(coordinates, 'min_points') && ...
        isequal(size(coordinates.min_points), [n_nodes, 3]), ...
        'min_points.mat must contain a 400-by-3 coordinate matrix.');
    if isfield(coordinates, 'label_set')
        assert(isequal(coordinates.label_set(:)', label_set), 'ROI coordinate labels must be 1-400.');
    end
    eu = euclidean_dis(coordinates.min_points);

    sc_file = fullfile(subject_dir, 'SC.mat'); % Weighted structural matrix, variable SC_L2R_order.
    binary_file = fullfile(subject_dir, 'SC_bi.mat'); % Binary structural matrix, variable SC_bi.
    pathlength_file = fullfile(subject_dir, 'pathlength.mat'); % Shortest-path matrix, variable sp.
    communicability_file = fullfile(subject_dir, 'communicability.mat'); % Communicability matrix, variable co.
    distance_file = fullfile(subject_dir, 'euclidean_dis.mat'); % Pairwise ROI distance matrix, variable eu.
    save(sc_file, 'SC_L2R_order', 'label_set');
    save(binary_file, 'SC_bi', 'label_set');
    save(pathlength_file, 'sp', 'label_set');
    save(communicability_file, 'co', 'label_set');
    save(distance_file, 'eu', 'label_set');
    fprintf('%s is done.\n', subject_id);
end
