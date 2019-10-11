function [MeanNegLogLike] = Minimize_This(CurrentParameters,Options,CDFData_Exp,PHIndex,FigureHandles,HexagonData)

disp(strcat('          ...Testing X = ',num2str(CurrentParameters)))
NumberPHValues = length(Options.pHValues);

% For the current set of parameters, evaluate the model at each pH value
for h = 1:NumberPHValues
    
    % If we are re-using the previous equilibration run, then set the output of
    % the previous simulation as the input of the current simulation 
    % Note: we reuse the equilibration run only for each set of pH
    % values, not across different parameters. So if you are
    % scanning parameters, then a separate equilibration run will
    % occur for each set of parameters. 
    SimInput.UsePriorEqRun = 'n';
    if strcmp(Options.ReuseEq,'y')
        if h ~= 1
            SimInput = SimOutput;
            SimInput.UsePriorEqRun = 'y';
        end
    end

    % Record the current pH value
    SimInput.pH = Options.pHValues(h);

    % Set up the rate constants and other parameters which are particular to the model you are using
    [SimInput] = Setup_Rate_Constants(CurrentParameters,SimInput,Options);

    % Initialize the simulation using the given set of parameters. Then
    % run the simulation. The equilibration run will also be performed
    % prior to the simulation run, if selected. 
        [SimOutput] = Initialize_And_Run_Model_GF(Options,CurrentParameters,...
            HexagonData,SimInput,FigureHandles);

        % If any fusion events were observed, calculate the CDF and randomness parameter
            if ~isempty(SimOutput.FusionWaitTimes)
                [CumX, CumY] = Calculate_CDF(SimOutput.FusionWaitTimes);
                    CumYNormalized = CumY/max(CumY);

                [RandomnessParameter,Nmin] = Calculate_Randomness_Parameter(SimOutput.FusionWaitTimes);
            else
                CumX = 0;
                CumY = 0;
                CumYNormalized = 0;
                RandomnessParameter = NaN;
                Nmin = NaN;
            end

        % Record the efficiency value of the first pH value in the
        % sequence, to which we will calculate the normalized efficiency 
        if h == 1
            FirstEfficiency = CumY(end)/SimOutput.NumberVirions;
        end

    % Compile data to save across pH values
%         ModelData(h).SimOutput = SimOutput;
        ModelData(h).CurrentParameters = CurrentParameters;
        ModelData(h).pH = SimOutput.pH;
        ModelData(h).CDFData.CumX = CumX;
        ModelData(h).CDFData.FusionWaitTimes = SimOutput.FusionWaitTimes;
        ModelData(h).CDFData.CumY = CumY;
        ModelData(h).CDFData.CumYNormalized = CumYNormalized;
        ModelData(h).CDFData.RandomnessParameter = RandomnessParameter;
        ModelData(h).CDFData.Nmin = Nmin;
        ModelData(h).CDFData.NumberVirions = SimOutput.NumberVirions;
        ModelData(h).CDFData.Efficiency = CumY(end)/SimOutput.NumberVirions;
        ModelData(h).CDFData.EfficiencyNorm = CumY(end)/FirstEfficiency;
%         ModelData(h).HelpfulInfo.PSetNumber = p;
        ModelData(h).HelpfulInfo.NumberMonomers = length(HexagonData.MonomerInfoLibrary);
        
end

% Calculate negative mean log likelihood of the model with the current
% parameters, given our experimental data
    [MeanNegLogLike] = Calculate_Likelihood_Global_Fit(ModelData,CDFData_Exp,PHIndex,Options,FigureHandles);
end