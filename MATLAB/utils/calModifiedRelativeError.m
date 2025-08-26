function relativeError = calModifiedRelativeError(valueTrue, valueIdentified, gamma)
if nargin < 3
    gamma = 0.01;
end
relativeError = abs(valueTrue-valueIdentified)/( gamma + abs(valueTrue));