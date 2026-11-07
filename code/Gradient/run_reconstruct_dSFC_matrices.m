% RUN_RECONSTRUCT_DSFC_MATRICES Reconstruct symmetric Schaefer 400 coupling matrices.
% Written by Junrui Zhang, 2025
%
% Input: all_matrixes.mat with all_matrixes [subjects x 79800].
% Output: 400ROI.mat containing matrices, a cell array of [400 x 400] matrices.
% Method: MATLAB column-major upper-triangle order; zero diagonal.
% Dependencies: MATLAB.
% Usage: Set data_root below; run run_compute_dSFC.m first.

clear;
clc;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
input_file = fullfile(data_root, 'results', 'dSFC', 'all_matrixes.mat'); % Subject-level upper-triangle coupling vectors.
output_dir = fullfile(data_root, 'results', 'Gradient'); % Output directory for reconstructed matrices and group gradients.
output_file = fullfile(output_dir, '400ROI.mat'); % Cell array of reconstructed Schaefer 400 coupling matrices.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the Gradient module.
n_nodes = 400;  % the number of parcels
label_set = 1:n_nodes;

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
input_data = load(input_file);
assert(isfield(input_data, 'all_matrixes'), 'The input must contain all_matrixes.');
matrices = reconstruct_symmetric_matrices(input_data.all_matrixes, n_nodes);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
if isfield(input_data, 'subject_ids')
    subject_ids = input_data.subject_ids;
    assert(numel(subject_ids) == numel(matrices), 'Subject metadata does not match the matrix count.');
    save(output_file, 'matrices', 'subject_ids', 'label_set', '-v7.3');
else
    save(output_file, 'matrices', 'label_set', '-v7.3');
end
fprintf('Reconstructed %d subject matrices.\n', numel(matrices));
