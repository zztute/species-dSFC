% RUN_COMPUTE_GROUP_GRADIENTS Fit diffusion-map gradients to group coupling.
% Written by Junrui Zhang, 2025
%
% Input: 400ROI.mat containing subject-specific [400 x 400] coupling matrices.
% Output: group_gradient.mat containing gm, explained_ratio, and mean_z_matrix.
% Method: Average Fisher-z matrices; normalized-angle kernel; diffusion mapping.
% Dependencies: BrainSpace MATLAB (GradientMaps).
% Usage: Set data_root below; run run_reconstruct_dSFC_matrices.m first.

clear;
clc;

%% Paths and parameters
data_root = '<PATH_TO_PROJECT_DATA>'; % Root directory containing the human example data.
gradient_dir = fullfile(data_root, 'results', 'Gradient'); % Directory containing reconstructed matrices and gradient outputs.
input_file = fullfile(gradient_dir, '400ROI.mat'); % Subject-level symmetric Schaefer 400 coupling matrices.
output_file = fullfile(gradient_dir, 'group_gradient.mat'); % Group GradientMaps object, eigenvalue proportions, and group mean Fisher-z matrix.
module_dir = fileparts(mfilename('fullpath')); % Directory containing this script.
addpath(fullfile(module_dir, 'function')); % Auxiliary functions for the Gradient module.
n_nodes = 400;  % the number of parcels
n_components = 10;
sparsity = 90;
label_set = 1:n_nodes;

assert(isfolder(data_root), 'Set data_root to an existing project-data directory.');
assert(exist('GradientMaps', 'class') == 8 || exist('GradientMaps', 'file') == 2, ...
    'Add BrainSpace MATLAB to the MATLAB path; see README.md.');
input_data = load(input_file, 'matrices');
assert(isfield(input_data, 'matrices') && iscell(input_data.matrices) && ...
    ~isempty(input_data.matrices), '400ROI.mat must contain a nonempty matrices cell array.');
matrices = input_data.matrices;
n_subjects = numel(matrices);
sum_matrix = zeros(n_nodes);

%% Average transformed subject matrices
for subject = 1:n_subjects
    matrix = matrices{subject};
    assert(isequal(size(matrix), [n_nodes, n_nodes]) && ...
        norm(matrix - matrix', 'fro') < 1e-8, ...
        'Each subject must provide a symmetric 400-by-400 matrix.');
    sum_matrix = sum_matrix + fisherz(matrix);
end
mean_z_matrix = sum_matrix / n_subjects;

%% Fit group gradients
gm = GradientMaps('kernel', 'na', 'approach', 'dm', 'n_components', n_components);
gm = gm.fit(mean_z_matrix, 'sparsity', sparsity);
lambda_vals = gm.lambda{1, 1};
explained_ratio = lambda_vals ./ sum(lambda_vals);
save(output_file, 'gm', 'explained_ratio', 'mean_z_matrix', 'label_set');
fprintf('Group gradient calculation complete.\n');
