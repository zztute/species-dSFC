% RECONSTRUCT_SYMMETRIC_MATRICES Restore matrices from upper-triangle vectors.
% Written by Junrui Zhang, 2025
%
% Inputs: edge_vectors [subjects x n_nodes*(n_nodes-1)/2]; scalar n_nodes.
% Output: matrices, a column cell array of symmetric matrices with zero diagonals.
% Ordering: MATLAB column-major indexing of triu(true(n_nodes), 1).
% Dependencies: MATLAB.

function matrices = reconstruct_symmetric_matrices(edge_vectors, n_nodes)
validateattributes(n_nodes, {'numeric'}, {'scalar', 'integer', '>=', 2, 'finite'});
validateattributes(edge_vectors, {'numeric'}, {'2d', 'nonempty', 'finite', 'real'});
mask = triu(true(n_nodes), 1);
assert(size(edge_vectors, 2) == nnz(mask), 'The edge-vector length does not match n_nodes.');
matrices = cell(size(edge_vectors, 1), 1);
for subject = 1:size(edge_vectors, 1)
    matrix = zeros(n_nodes);
    matrix(mask) = edge_vectors(subject, :);
    matrices{subject} = matrix + matrix';
end
end
