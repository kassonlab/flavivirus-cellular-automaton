function [] = Start_Fit_CA(ModelToUse,CurrentInstance,NumberInstances)

% - - - - - - - - - - - - - - - - - - - - -
% Input:
% Start_Fit_CA(ModelToUse,CurrentInstance,NumberInstances,SavePath)
%     ModelToUse = The name of the model which will be run, see
%     Setup_Rate_Constants.m. Current options are OffPath,
%     Linear2pHRateContants, Linear1pHRateContant  
% 
%     CurrentInstance = The number of the current instance of Matlab if you
%     are running multiple instances as a crude parallelization. 
% 
%     NumberInstances = The total number of instances of Matlab you are
%     running if performing a crude parallelization. 

% Output:
% Solution displayed at command prompt

% By Bob Rawle, Kasson Lab, University of Virginia, 2017
% - - - - - - - - - - - - - - - - - - - - -

close all
tic
% Set up options, including the initial guesses for the fit
    [Options] = Setup_Options_Fit_CA();
    
% load the experimental data
    ExpData = open(Options.ExpDataFilename);
    CDFData_Exp = ExpData.CDFData;
    
    % Create pH index for easy reference, also use relative efficiencies if needed
    for e = 1:length(CDFData_Exp)
        PHIndex(e) = CDFData_Exp(e).pH;
        if strcmp(Options.UseRelEfficiencies,'y')
            CDFData_Exp(e).EfficiencyCorrected = CDFData_Exp(e).EfficiencyBefore/Options.EffMax;
        else
            CDFData_Exp(e).EfficiencyCorrected = CDFData_Exp(e).EfficiencyBefore;
        end
    end

% Set up figure handles
    if strcmp(Options.DisplayFigures,'y')
        [FigureHandles] = Setup_Figures_Fit_CA(Options);
    else
        FigureHandles = [];
    end

% Use preset hexagon file, or set up the hexagon, defining dimers as well as possible trimers for each monomer
    if ~isempty(Options.PresetHexagonFilename)
        PresetHexagonData = load(Options.PresetHexagonFilename);
        HexagonData = PresetHexagonData.HexagonData;
    else
        [HexagonData] = Setup_Hexagon();
    end

% Record which model we are using
        Options.ModelToUse = ModelToUse;
    
% Split up the initial guesses if you are doing a crude parallelization by
% opening multiple instances of Matlab
    NumberInitialGuesses = size(Options.GuessArray,1);
        if strcmp(Options.CrudeParallel,'y')
            NumGuessesPerInstance = floor(NumberInitialGuesses/NumberInstances);
            if CurrentInstance==NumberInstances
                GuessRange = (CurrentInstance-1)*NumGuessesPerInstance+1:NumberInitialGuesses;
            else
                GuessRange = (CurrentInstance-1)*NumGuessesPerInstance+1:CurrentInstance*NumGuessesPerInstance;
            end
        else 
            GuessRange = 1:NumberInitialGuesses;
        end
    
% For each set of initial guesses, run the fit to the experimental data
    for p = GuessRange
            CurrentGuess = Options.GuessArray(p,:);
            
            disp(' ');disp(' ');
            disp(strcat('----InitialGuess =',num2str(CurrentGuess),' -----'));
            
            % Perform the fit using the chosen algorithm
            if strcmp(Options.FitMethod,'Constrained-PatternSearch')
                LowerBound = Options.LowerBounds(p,:);
                UpperBound = Options.UpperBounds(p,:);
                disp(strcat('  LowerBound =',num2str(LowerBound)));
                disp(strcat('  UpperBound =',num2str(UpperBound)));
                
                    OptimizationOptions = psoptimset('TolX', 1e-5,'TolFun', 1e-6,'Display','iter',...
                        'PlotFcns', {@displaybestx,@psplotbestf},'MaxIter',200,...
                        'InitialMeshSize',.3*min(CurrentGuess),...
                        'ScaleMesh','Off','PollingOrder','Success','Cache','on');
                
                Min_Fun = @(CurrentGuess)Minimize_This(CurrentGuess,Options,CDFData_Exp,PHIndex,FigureHandles,HexagonData);

                [Solution,MeanNegLogLike] = patternsearch(Min_Fun,CurrentGuess,[],[],[],[],...
                    LowerBound,UpperBound,[],OptimizationOptions);
            end
            
        
            disp(strcat('Solution =',num2str(Solution)));
            disp(strcat('MeanNegLogLikelihood=',num2str(MeanNegLogLike)));
            
            % Evaluate and plot the fit solution
            Evaluate_And_Plot_Fit_Solution(Solution,Options,FigureHandles,CDFData_Exp,PHIndex,HexagonData);
        
    end
    
    toc
    disp('Thank You, Come Again')

end