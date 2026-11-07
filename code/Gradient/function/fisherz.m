% FISHERZ Apply the Fisher r-to-z transform element by element.
% Written by Junrui Zhang, 2025
%
% Input: R, a real finite correlation array with values strictly inside (-1, 1).
% Output: Z = atanh(R), preserving input dimensions.
% Behavior: Perfect correlations are rejected rather than clipped.
% Dependencies: MATLAB.

function Z = fisherz(R)
validateattributes(R, {'numeric'}, {'nonempty', 'finite', 'real'});
assert(all(abs(R(:)) < 1), ...
    'Fisher-z requires correlations strictly inside (-1, 1); inspect perfect correlations.');
Z = atanh(double(R));
end
