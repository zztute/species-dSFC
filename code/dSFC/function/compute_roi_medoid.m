% COMPUTE_ROI_MEDOID Select the observed point with minimum total distance.
% Written by Junrui Zhang, 2025
%
% Input: coordinates, a finite nonempty [voxels x 3] matrix.
% Output: medoid, the coordinate minimizing summed Euclidean distance.
% Method: Exact exhaustive search with linear temporary-memory usage.
% Dependencies: MATLAB.

function medoid = compute_roi_medoid(coordinates)
validateattributes(coordinates, {'numeric'}, {'2d', 'ncols', 3, 'nonempty', 'finite', 'real'});
coordinates = double(coordinates);
n_points = size(coordinates, 1);
total_distances = zeros(n_points, 1);
for point = 1:n_points
    differences = coordinates - coordinates(point, :);
    total_distances(point) = sum(sqrt(sum(differences.^2, 2)));
end
[~, minimum_index] = min(total_distances);
medoid = coordinates(minimum_index, :);
end
