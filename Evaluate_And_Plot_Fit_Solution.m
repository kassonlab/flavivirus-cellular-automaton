function [] = Evaluate_And_Plot_Fit_Solution(CurrentParams,Options,FigureHandles,CDFData_Exp,PHIndex,HexagonData)

    NumberPHValues = length(Options.pHValues);
    PlotCounter = 0;

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
    [SimInput] = Setup_Rate_Constants(CurrentParams,SimInput,Options);

    % Initialize the simulation using the given set of parameters. Then
    % run the simulation. The equilibration run will also be performed
    % prior to the simulation run, if selected. 
        [SimOutput] = Initialize_And_Run_Model_GF(Options,CurrentParams,...
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

        % Display the figures        
            if strcmp(Options.DisplayFigures,'y')
                % Plot the solution
                set(0,'CurrentFigure',FigureHandles.ModelPlot)
                hold on
                
                if strcmp(Options.NormalizeModelEfficiencies,'y')
                    Efficiency_Model(h) = CumY(end)/FirstEfficiency/SimOutput.NumberVirions;
                    plot(CumX,CumY/FirstEfficiency,'x')
                else
                    Efficiency_Model(h) = CumY(end)/SimOutput.NumberVirions;
                    plot(CumX,CumY/SimOutput.NumberVirions,'x')
                end

                set(0,'CurrentFigure',FigureHandles.ModelNormalizedPlot)
                hold on
                plot(CumX,CumYNormalized,'x')
                drawnow
                
                PlotCounter = PlotCounter + 1;
                LegendInfo_Model{1,h} = strcat('pH=',...
                    num2str(SimOutput.pH),'; Nmin=',num2str(Nmin,'%.2f'),'; Eff=',num2str(CumY(end)/SimOutput.NumberVirions,...
                        '%.2f'));

                % load up the experimental data for the current pH value
                ExpIndex = find(PHIndex ==SimOutput.pH);
                CumX_Exp = CDFData_Exp(ExpIndex).CumX;
                CumYNormalized_Exp = CDFData_Exp(ExpIndex).CumYNorm;
                Efficiency_Exp(h) = CDFData_Exp(ExpIndex).EfficiencyCorrected;
                CumY_Exp = CumYNormalized_Exp*Efficiency_Exp(h);

                % plot experimental data from the current pH value
                set(0,'CurrentFigure',FigureHandles.ExpPlot)
                plot(CumX_Exp,CumY_Exp,'o')
                hold on

                set(0,'CurrentFigure',FigureHandles.ExpNormalizedPlot)
                plot(CumX_Exp,CumYNormalized_Exp,'o')
                hold on

                PlotCounter = PlotCounter + 1;
                LegendInfo_Exp{1,h} = strcat('pH=',num2str(SimOutput.pH),'; Eff=',num2str(Efficiency_Exp(h),'%.2f'));

            end
    end


    
    if strcmp(Options.DisplayFigures,'y')
        % Make efficiency plot
            set(0,'CurrentFigure',FigureHandles.EfficiencyPlot)
            plot(Options.pHValues,Efficiency_Exp,'o-');
            hold on
            plot(Options.pHValues,Efficiency_Model,'o-');
            ylim([0 1])
            LegendInfo_Efficiency{1,1} = strcat('Exp Data');
            LegendInfo_Efficiency{1,2} = strcat('Model');
            legend(LegendInfo_Efficiency,'Location','best');
            xlabel('pH');
            ylabel('Efficiency');

        % Put legends on plots
        set(0,'CurrentFigure',FigureHandles.ModelPlot)
        legend(LegendInfo_Model,'Location','best');
        xlabel('Time');
        ylabel('Percent Fusion');
        title(strcat('Solution=',num2str(CurrentParams)))

        set(0,'CurrentFigure',FigureHandles.ModelNormalizedPlot)
        legend(LegendInfo_Model,'Location','best');
        xlabel('Time');
        ylabel('Normalized Fusion');
        title(strcat('Solution=',num2str(CurrentParams)))
        
        set(0,'CurrentFigure',FigureHandles.ExpPlot)
        legend(LegendInfo_Exp,'Location','best');
        xlabel('Time');
        ylabel('Percent Fusion');
        title(strcat('ExperimentalData'))

        set(0,'CurrentFigure',FigureHandles.ExpNormalizedPlot)
        legend(LegendInfo_Exp,'Location','best');
        xlabel('Time');
        ylabel('Normalized Fusion');
        title(strcat('ExperimentalData'))
    end

end