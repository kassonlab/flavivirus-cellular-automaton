function [ProbDensityFused] = Estimate_Kernal_Density(FusionWaitTimes_Model,ExpWaitTimes,Options)

% Implement boundary reflection of data to better estimate the density at the edges of the data set
    TimeCutoff = Options.TotalTime;
    ReflectionWindow = 100;
        % Time in seconds over which we will reflect the data across the boundary

    ReflectionValuesNeg = FusionWaitTimes_Model(FusionWaitTimes_Model < ReflectionWindow) * -1;
    ReflectionValuesPos = FusionWaitTimes_Model(FusionWaitTimes_Model > TimeCutoff - ReflectionWindow);
    ReflectionToAdd = TimeCutoff - ReflectionValuesPos;
    ReflectionValuesPos = ReflectionValuesPos + 2*ReflectionToAdd;
    
    FusionWaitTimes_Model = [ReflectionValuesNeg FusionWaitTimes_Model ReflectionValuesPos];

% Calculate the kernel density estimate and then re-normalize across the original support (not including the reflected data)
[ProbDensityFused, ~, Bandwidth] = ksdensity(FusionWaitTimes_Model,ExpWaitTimes,'Support','unbounded',...
    'Kernel','normal');

ProbDensityFused= ProbDensityFused/trapz(ExpWaitTimes,ProbDensityFused);

end