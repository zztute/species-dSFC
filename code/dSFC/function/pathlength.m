% PATHLENGTH Calculate shortest-path lengths in an unweighted network.
% Written by Junrui Zhang, 2025
%
% Input: SC_bi, a binary square adjacency matrix.
% Output: sp, a shortest-path matrix; disconnected pairs have Inf entries.
% Dependencies: Brain Connectivity Toolbox (distance_bin).

function sp = pathlength(SC_bi)
validateattributes(SC_bi, {'numeric', 'logical'}, {'2d', 'square', 'nonempty', 'finite', 'real'});
assert(all(SC_bi(:) == 0 | SC_bi(:) == 1), 'SC_bi must be a binary adjacency matrix.');
assert(exist('distance_bin', 'file') == 2, 'Add Brain Connectivity Toolbox to the MATLAB path.');
sp = distance_bin(double(SC_bi));
end
