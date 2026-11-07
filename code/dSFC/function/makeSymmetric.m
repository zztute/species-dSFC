% MAKESYMMETRIC Convert FC matrix cells from sparse to full storage.
% Written by Junrui Zhang, 2025
%
% Input: B, a cell array of FC matrices.
% Output: B_full, a column cell array containing the same matrix values.
% Behavior: The supplied function name is retained; no symmetrization is applied.
% Dependencies: MATLAB.

function B_full = makeSymmetric(B)
assert(iscell(B) && ~isempty(B), 'FCM.Matrix must be a nonempty cell array.');
n_cells = numel(B);
B_full = cell(n_cells, 1);
for index = 1:n_cells
    validateattributes(B{index}, {'numeric'}, {'2d', 'square', 'nonempty', 'real'});
    B_full{index} = full(B{index});
end
end
