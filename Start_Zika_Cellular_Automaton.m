function [] = Start_Zika_Cellular_Automaton(ModelToUse,CurrentInstance,NumberInstances,SavePath)

% - - - - - - - - - - - - - - - - - - - - -
% Input:
% Start_Zika_Cellular_Automaton(ModelToUse,CurrentInstance,NumberInstances,SavePath)
%     ModelToUse = The name of the model which will be run, see
%     Setup_Rate_Constants.m. Current options are OffPath,
%     Linear2pHRateContants, Linear1pHRateContant  
% 
%     CurrentInstance = The number of the current instance of Matlab if you
%     are running multiple instances as a crude parallelization. 
% 
%     NumberInstances = The total number of instances of Matlab you are
%     running if performing a crude parallelization. 
% 
%     SavePath = Directory where the output .mat file will be saved

% Output:
% If chosen, the cellular automaton data will be saved as a .MAT file.

% By Bob Rawle, Kasson Lab, University of Virginia, 2017, based on the
% algorithm described in Chao et al., 2014, eLife, doi:10.7554/eLife.04389
% - - - - - - - - - - - - - - - - - - - - -
tic
% Set up options, including the parameters to scan
    [Options] = Setup_Options();

    if strcmp(Options.SaveData,'y') 
        if ~isempty(SavePath)
            SaveFolderDir = SavePath;
        else
            SaveFolderDir = uigetdir('Choose the directory where data folder will be saved');
        end
    else
        SaveFolderDir = [];
    end
    
    % Set up larger parameter scan, if chosen
    if ischar(Options.ParameterArray) && strcmp(Options.ParameterArray,'Large')
        Options = Setup_Large_Parameter_Scan(Options,SaveFolderDir);
    end
        
        if isnan(Options.ParameterArray)
            NumberParameterCycles = 1;
        else 
            NumberParameterCycles = size(Options.ParameterArray,1);
        end
        
        % Define parameter range, which may change if you are running
        % multiple instances of Matlab as a crude parallelization 
        if strcmp(Options.Parallel,'y')
            NumCyclesPerInstance = floor(NumberParameterCycles/NumberInstances);
            if CurrentInstance==NumberInstances
                ParameterRange = (CurrentInstance-1)*NumCyclesPerInstance+1:NumberParameterCycles;
            else
                ParameterRange = (CurrentInstance-1)*NumCyclesPerInstance+1:CurrentInstance*NumCyclesPerInstance;
            end
        else 
            ParameterRange = 1:NumberParameterCycles;
        end
        
% Record which model we are using
    Options.ModelToUse = ModelToUse;


% Use preset hexagon file, or set up the hexagon, defining dimers as well as possible trimers for each monomer
    if ~isempty(Options.PresetHexagonFilename)
        PresetHexagonData = load(Options.PresetHexagonFilename);
        HexagonData = PresetHexagonData.HexagonData;
    else
        [HexagonData] = Setup_Hexagon();
    end

% Set up figure handles
    if strcmp(Options.DisplayFigures,'y')
        [FigureHandles] = Setup_Figures(Options);
    else
        FigureHandles = [];
    end

% Initialize other variables
    CycleCounter = 0;
    NumberPHValues = length(Options.pHValues);
    
% For each set of parameters to scan and pH values, run the simulation and compile the data
    for p = ParameterRange
        % Reset the data to save for each parameter set
            DataToSave = [];
            
        for h = 1:NumberPHValues
            CurrentParameters = Options.ParameterArray(p,:);
            CycleCounter = (p-1)*NumberPHValues+h;
            
            % If we are re-using the previous equilibration run, then set the output of
            % the previous simulation as the input of the current simulation 
            % Note: we reuse the equilibration run only for each set of pH
            % values, not across different parameters. So if you are
            % scanning parameters, then a separate equilibration run will
            % occur for each set of parameters.   
            if strcmp(Options.ReuseEq,'y')
                if h == 1
                    SimInput.UsePriorEqRun = 'n';
                else
                    SimInput = SimOutput;
                    SimInput.UsePriorEqRun = 'y';
                end
            else 
                SimInput.UsePriorEqRun = 'n';
            end
            
            SimInput.pH = Options.pHValues(h);


            disp(' ');disp(' ');
            disp(strcat(' *** Simulation_',num2str(CycleCounter),'/',num2str(NumberParameterCycles*NumberPHValues),' ***'));                
            disp(strcat('----Params =',num2str(CurrentParameters),'; pH = ',num2str(SimInput.pH),' -----'));

            % Set up the rate constants and other parameters which are particular to the model you are using
                [SimInput] = Setup_Rate_Constants(CurrentParameters,SimInput,Options);


            % Initialize the simulation using the given set of parameters. Then
            % run the simulation.
                [SimOutput] = Initialize_And_Run_Model(Options,CurrentParameters,...
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

            % Plot the data
            if strcmp(Options.DisplayFigures,'y')
                set(0,'CurrentFigure',FigureHandles.MainPlot)
                hold on
                plot(CumX,CumY/SimOutput.NumberVirions*100,'o')

                set(0,'CurrentFigure',FigureHandles.NormalizedPlot)
                hold on
                plot(CumX,CumYNormalized,'o')
                drawnow

                LegendInfo{1,CycleCounter} = strcat('Param=',num2str(CurrentParameters),'; pH=',...
                    num2str(SimInput.pH),'; Nmin=',num2str(Nmin,'%.2f'));
            end

            % Compile data to save across pH values
            if strcmp(Options.SaveData,'y')
                DataToSave(h).SimOutput = SimOutput;
                DataToSave(h).CurrentParameters = CurrentParameters;
                DataToSave(h).CDFData.CumX = CumX;
                DataToSave(h).CDFData.CumY = CumY;
                DataToSave(h).CDFData.CumYNormalized = CumYNormalized;
                DataToSave(h).CDFData.RandomnessParameter = RandomnessParameter;
                DataToSave(h).CDFData.Nmin = Nmin;
                DataToSave(h).CDFData.NumberVirions = SimOutput.NumberVirions;
                DataToSave(h).CDFData.Efficiency = CumY(end)/SimOutput.NumberVirions*100;
                DataToSave(h).HelpfulInfo.PSetNumber = p;
                DataToSave(h).HelpfulInfo.NumberMonomers = length(HexagonData.MonomerInfoLibrary);
            end

        end
        
        % Save data for each parameter set to a separate file
        if strcmp(Options.SaveData,'y')
            save(strcat(SaveFolderDir,'/',Options.FileLabel,'_Params=',num2str(CurrentParameters),'.mat'),'DataToSave')
        end
    end

    % Put legends on plots
    if strcmp(Options.DisplayFigures,'y')
        set(0,'CurrentFigure',FigureHandles.MainPlot)
        legend(LegendInfo,'Location','best');
        xlabel('Time');
        ylabel('Percent Fusion');

        set(0,'CurrentFigure',FigureHandles.NormalizedPlot)
        legend(LegendInfo,'Location','best');
        xlabel('Time');
        ylabel('Normalized Fusion');
    end

    toc
    disp('Thank You, Come Again')

end