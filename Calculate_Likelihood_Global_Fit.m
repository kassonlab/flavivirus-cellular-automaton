function [MeanNegLogLike] = Calculate_Likelihood_Global_Fit(ModelData,CDFData_Exp,PHIndex,Options,FigureHandles)
NumberpHValues = length(ModelData);

for h = 1:NumberpHValues
    % load up the data from the model for the current pH value
    pH_Model = ModelData(h).pH;
    CumX_Model = ModelData(h).CDFData.CumX;
    FusionWaitTimes_Model = ModelData(h).CDFData.FusionWaitTimes;
    CumYNormalized_Model = ModelData(h).CDFData.CumYNormalized;
    CumY_Model = ModelData(h).CDFData.CumY;
    
    if strcmp(Options.NormalizeModelEfficiencies,'y')
        Efficiency_Model = ModelData(h).CDFData.EfficiencyNorm;
    else 
        Efficiency_Model = ModelData(h).CDFData.Efficiency;
    end    
    
        % Record the gating efficiency for the gate test later on
        if strcmp(Options.LikelihoodMethod,'Kinetics,GatedEfficiencies')  &&...
                pH_Model == Options.GatePHValues(1)
            GatingEfficiencies(1) = Efficiency_Model;
        elseif strcmp(Options.LikelihoodMethod,'Kinetics,GatedEfficiencies')  &&...
                pH_Model == Options.GatePHValues(2)
            GatingEfficiencies(2) = Efficiency_Model;
        end

        if strcmp(Options.Diagnostics,'y')
            set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
            plot(CumX_Model,PDF_Model)
        end
    
    % load up the experimental data for the current pH value
    ExpIndex = find(PHIndex ==pH_Model);
    ExpWaitTimes = CDFData_Exp(ExpIndex).SortedpHtoFList;
    ExpEfficiency = CDFData_Exp(ExpIndex).EfficiencyCorrected;
    ExpNumberFused = length(ExpWaitTimes);
    ExpNumberNotFused = round(ExpNumberFused/ExpEfficiency) - ExpNumberFused;
        if ExpNumberNotFused < 0
            ExpNumberNotFused = 0;
        end
    ExpNumberTotal = ExpNumberFused + ExpNumberNotFused;
    
    % calculate the mean log likelihood for the current pH value. The
    % probability density for each experimental data point of a virus that
    % fused is calculated using the kernel density estimate from the model data
    if ~isempty(FusionWaitTimes_Model)
        [ProbDensityFused] = Estimate_Kernal_Density(FusionWaitTimes_Model,ExpWaitTimes,Options);
    else
        ProbDensityFused = 1;
            % If no fusion events were observed in the model, set the
            % probability density to 1, which will result in a 0 in the log
            % likelihood function below.  
    end
        
    if strcmp(Options.LikelihoodMethod,'KineticsAndEfficiencies')
        if Efficiency_Model == 1
            LogLikePH = sum(Options.KineticWeightConstant*log(ProbDensityFused)) + ExpNumberFused*log(Efficiency_Model);
        else
            LogLikePH = sum(Options.KineticWeightConstant*log(ProbDensityFused)) + ExpNumberFused*log(Efficiency_Model) + ExpNumberNotFused*log(1-Efficiency_Model);
            % Note that as Efficiency_Model goes to zero, this will explode to negative infinity
        end
        MeanLogLikePH(h) = LogLikePH/ExpNumberTotal;
    elseif strcmp(Options.LikelihoodMethod,'Kinetics,GatedEfficiencies') 
        LogLikePH = sum(log(ProbDensityFused));
        MeanLogLikePH(h) = LogLikePH/ExpNumberFused;
    end
    
end

    % calculate the sum of all mean log likelihoods for each pH value. Multiply
    % by -1 to yield the negative mean log likelihood for the entire data set.
    MeanNegLogLike = -1*sum(MeanLogLikePH);
    
    %Gate test
    if strcmp(Options.LikelihoodMethod,'Kinetics,GatedEfficiencies') 
        if GatingEfficiencies(1) - GatingEfficiencies(2) >= Options.GateValue ...
                && GatingEfficiencies(1) >= Options.HighEffValueCutoff
        else
            MeanNegLogLike = NaN;
        end
    end
end