function [factor_downsample,tsampIdxPeriodic_downsample, ...
    numTPeriodic_downsample,deltaT_downsample] = getDownsample( ...
    numT,deltaT,factor_downsample)

% NOTE: The last timepoint is actually the start of the period, so exclude
tsampIdxPeriodic_downsample = 1:factor_downsample:(numT-1);
numTPeriodic_downsample = numel(tsampIdxPeriodic_downsample);
deltaT_downsample = factor_downsample*deltaT;
end