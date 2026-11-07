% FCN_COMMUNICABILITY Calculate unweighted network communicability.
% Written by Junrui Zhang, 2025
%
% Input: A, a binary square adjacency matrix.
% Output: G = expm(A), the matrix exponential of the adjacency matrix.
% Dependencies: MATLAB.
% Original function author: Bratislav Misic (Bratislav Mišić).
% Adaptation: Header and function name organized for the supplied workflow.
% Reference: Estrada, E. and Hatano, N. (2008). Communicability in complex
% networks. Physical Review E, 77, 036111. doi:10.1103/PhysRevE.77.036111.

function G = fcn_communicability(A)
validateattributes(A, {'numeric', 'logical'}, {'2d', 'square', 'nonempty', 'finite', 'real'});
assert(all(A(:) == 0 | A(:) == 1), 'A must be an unweighted binary adjacency matrix.');
G = expm(double(A));
end
