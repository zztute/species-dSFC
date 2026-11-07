% RUN_COMPUTE_DSFC Calculate dynamic structure-function coupling.
% Written by Junrui Zhang, 2025
%
% Inputs: DynamicBC FCM structures and three structural predictors per subject.
% Outputs: all_matrixes [subjects x 79800] and all_vsfc [subjects x 400].
% Method: Adjusted R-squared of regional linear models in each FC window.
% Dependencies: Statistics and Machine Learning Toolbox.
% Usage: Set data_root below; calculate structural predictors and DynamicBC FC first.
% Provenance: Supplied version optimized by Delmor, based on Zhang Zhe (2023).

clear;
clc;
tic;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
fc_root = fullfile(data_root, 'DynamicFC'); % Subject folders, each containing one DynamicBC MAT file with FCM.
sc_root = fullfile(data_root, 'SC', 'Subj'); % Subject folders containing pathlength.mat, communicability.mat, and euclidean_dis.mat.
output_dir = fullfile(data_root, 'results', 'dSFC'); % Output directory for dynamic structure-function coupling summaries.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the dSFC module.
n_nodes = 400; % the number of parcels
label_set = 1:n_nodes;
upper_mask = triu(true(n_nodes), 1);

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
assert(exist('fitlm', 'file') == 2, 'Install Statistics and Machine Learning Toolbox.');
fc_subject_ids = get_subject_ids(fc_root);
sc_subject_ids = get_subject_ids(sc_root);
subject_ids = intersect(fc_subject_ids, sc_subject_ids);
assert(~isempty(subject_ids), 'No matching structural and dynamic-FC subject folders were found.');
if numel(subject_ids) < numel(fc_subject_ids) || numel(subject_ids) < numel(sc_subject_ids)
    warning('Only subjects present in both input directories are analyzed.');
end
fprintf('Found %d common subjects.\n', numel(subject_ids));
all_matrixes = zeros(numel(subject_ids), nnz(upper_mask));
all_vsfc = zeros(numel(subject_ids), n_nodes);

%% Fit regional models in each dynamic-FC window
for subject = 1:numel(subject_ids)
    subject_id = subject_ids{subject};
    fprintf('%s is processing...\n', subject_id);
    fc_subject_dir = fullfile(fc_root, subject_id); % Directory containing this subject's DynamicBC FC output.
    fc_files = dir(fullfile(fc_subject_dir, '*.mat')); % Candidate DynamicBC MAT files; exactly one is required per subject.
    assert(isscalar(fc_files), ...
        'Subject %s must have exactly one DynamicBC MAT file in its FC folder.', subject_id);
    fc_file = fullfile(fc_subject_dir, fc_files(1).name); % DynamicBC output MAT file containing FCM.Matrix and FCM.time_alignment.
    fc_data = load(fc_file, 'FCM');
    assert(isfield(fc_data, 'FCM') && isstruct(fc_data.FCM) && ...
        isfield(fc_data.FCM, 'Matrix') && isfield(fc_data.FCM, 'time_alignment'), ...
        'The DynamicBC input must contain FCM.Matrix and FCM.time_alignment.');
    FCM = fc_data.FCM;
    dFC = makeSymmetric(FCM.Matrix);
    n_windows = numel(FCM.time_alignment);
    assert(n_windows == numel(dFC) && n_windows >= 2, ...
        'DynamicBC FC matrices and window-alignment entries must have matching lengths.');

    sc_subject_dir = fullfile(sc_root, subject_id); % Structural predictors for the matching subject.
    pathlength_file = fullfile(sc_subject_dir, 'pathlength.mat'); % Shortest-path matrix MAT file, variable sp.
    communicability_file = fullfile(sc_subject_dir, 'communicability.mat'); % Communicability matrix MAT file, variable co.
    distance_file = fullfile(sc_subject_dir, 'euclidean_dis.mat'); % ROI Euclidean-distance matrix MAT file, variable eu.
    pathlength_data = load(pathlength_file, 'sp');
    communicability_data = load(communicability_file, 'co');
    distance_data = load(distance_file, 'eu');
    assert(isfield(pathlength_data, 'sp') && isfield(communicability_data, 'co') && ...
        isfield(distance_data, 'eu'), 'Structural predictor variable names must be sp, co, and eu.');
    predictors = {pathlength_data.sp, communicability_data.co, distance_data.eu};
    for predictor = 1:numel(predictors)
        assert(isequal(size(predictors{predictor}), [n_nodes, n_nodes]) && ...
            all(isfinite(predictors{predictor}(:))), ...
            'Structural predictors must be finite 400-by-400 matrices; check disconnected SC nodes.');
    end

    dsfc = zeros(n_windows, n_nodes);
    for window = 1:n_windows
        window_fc = dFC{window};
        assert(isequal(size(window_fc), [n_nodes, n_nodes]) && all(isfinite(window_fc(:))), ...
            'Each DynamicBC FC window must contain a finite 400-by-400 matrix.');
        assert(norm(window_fc - window_fc', 'fro') <= 1e-8 * max(1, norm(window_fc, 'fro')), ...
            'DynamicBC FC matrices must already be symmetric; provide both triangles.');
        for region = 1:n_nodes
            y = window_fc(:, region);
            x1 = pathlength_data.sp(:, region);
            x2 = communicability_data.co(:, region);
            x3 = distance_data.eu(:, region);
            x = zscore([x1, x2, x3]);
            lm = fitlm(x, y, 'Exclude', region);
            dsfc(window, region) = lm.Rsquared.Adjusted;
        end
    end
    coupling_matrix = corr(dsfc);
    assert(all(isfinite(coupling_matrix(:))), ...
        'Undefined coupling correlations for subject %s; inspect regional model scores.', subject_id);
    all_matrixes(subject, :) = coupling_matrix(upper_mask)';
    all_vsfc(subject, :) = var(dsfc, 0, 1);
    fprintf('%s is done.\n', subject_id);
end

%% Save subject-level summaries
if ~isfolder(output_dir)
    mkdir(output_dir);
end
matrices_file = fullfile(output_dir, 'all_matrixes.mat'); % Upper-triangle coupling vectors, subject IDs, and ROI labels.
variability_file = fullfile(output_dir, 'all_vsfc.mat'); % Regional coupling-score variance, subject IDs, and ROI labels.
save(matrices_file, 'all_matrixes', 'subject_ids', 'label_set', '-v7.3');
save(variability_file, 'all_vsfc', 'subject_ids', 'label_set', '-v7.3');
toc;
