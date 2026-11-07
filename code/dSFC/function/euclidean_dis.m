% EUCLIDEAN_DIS Calculate pairwise distances between ROI coordinates.
% Written by Junrui Zhang, 2025
%
% Input: min_points, a finite [ROIs x 3] coordinate matrix.
% Output: eu, a symmetric [ROIs x ROIs] matrix with a zero diagonal.
% Units: Same coordinate units as min_points.
% Dependencies: Statistics and Machine Learning Toolbox (pdist, squareform).

function eu = euclidean_dis(min_points)
validateattributes(min_points, {'numeric'}, {'2d', 'ncols', 3, 'nonempty', 'finite', 'real'});
eu = squareform(pdist(double(min_points)));
end
