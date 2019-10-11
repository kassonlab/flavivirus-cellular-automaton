function [RandomnessParameter,Nmin] = Calculate_Randomness_Parameter(WaitTimeList)

    MeanWaitingTimeSquared = (mean(WaitTimeList))^2;
    MeanOfSquared = mean(WaitTimeList.^2);

    RandomnessParameter = (MeanOfSquared - MeanWaitingTimeSquared)/MeanWaitingTimeSquared;
    Nmin = 1/RandomnessParameter;

end