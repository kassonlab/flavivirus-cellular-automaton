function [CumProbLookup,ProbLookup,TimeStep] = Det_Time_Step(k,CoopFactor,Options,EQRun)
% Determine the time step which will give us the specified maximum probability

if strcmp(EQRun,'y')
    ProbCutoff = Options.ProbCutoff_Eq;
else
    ProbCutoff = Options.ProbCutoff;
end

% Make an initial guess based on a single step CDF of the maximum rate constant (non-diagonal)
    k_guess = k;
    k_guess(1,2) = CoopFactor*k_guess(1,2);
    Index = eye(size(k_guess));
    k_max = max(k_guess(~Index));
    InitialGuess = -log(-(ProbCutoff - 1))/k_max;
    
% Determine the time step
    UpperBound = InitialGuess*100;
    LowerBound = InitialGuess/100;
    OptimizationOptions = optimset('TolX', 1e-16,'TolFun', 1e-14,'Algorithm', 'interior-point','Display','none');

    TimeStep = fmincon(@Time_Step_Min_Fun,InitialGuess,[],[],[],[],...
        LowerBound,UpperBound,[],OptimizationOptions,k,CoopFactor,Options,ProbCutoff);

    if TimeStep < Options.MinAdaptiveTimeStep
        TimeStep = Options.MinAdaptiveTimeStep;
    end
    
% Now calculate the probability lookup table with the determined time step
    [CumProbLookup,ProbLookup] = Compute_Probability_Lookup_Table(k,TimeStep,CoopFactor,Options);
    
    %if strcmp(Options.Diagnostics,'y')
        OverallMaxProb = Find_Max_Prob(ProbLookup);
        disp(strcat('   Adaptive Time Step = ',num2str(TimeStep),' s'))
        disp(strcat('   Max Probability= ',num2str(OverallMaxProb)))
    %end
    
end

function Error = Time_Step_Min_Fun(TimeStep,k,CoopFactor,Options,ProbCutoff)

[~,ProbLookup] = Compute_Probability_Lookup_Table(k,TimeStep,CoopFactor,Options);

OverallMaxProb = Find_Max_Prob(ProbLookup);

Error = abs(OverallMaxProb - ProbCutoff);
end

function OverallMaxProb = Find_Max_Prob(ProbLookup)
    OverallMaxProb = 0;
    
    % Find all possible combinations of geometrical states
    ProbNumDim = ndims(ProbLookup);
    PossibleValues = cell(1,ProbNumDim);
    for i =  1:ProbNumDim
        PossibleValues{1,i} = [ 1:size(ProbLookup,i)];
    end

    AllProbStateCombos = allcomb(PossibleValues{1,1:ProbNumDim});
    
    % Examine the probability matrices for all geometrical state
    % possibilities and determine the largest non-diagonal probability 
    for j = 1:size(AllProbStateCombos,1)
        CurrProbMatrix = ProbLookup{AllProbStateCombos(j,:)};
        Index = eye(size(CurrProbMatrix));
        MaxProb = max(CurrProbMatrix(~Index));
        if MaxProb > OverallMaxProb
            OverallMaxProb = MaxProb;
        end
    end
end