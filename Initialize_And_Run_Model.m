function [SimOutput] = Initialize_And_Run_Model(Options,CurrentParameters,...
    HexagonData,SimInput,FigureHandles)

% Open up Options, etc. for easier reference
    NumberVirions = Options.NumberVirions;
    NumberMonomers = length(HexagonData.MonomerInfoLibrary);

% Set up the E state matrix (contains the current state of all monomers in
% all viruses).By default, we say that all start in state 1
    EStates = ones(NumberVirions,NumberMonomers);
    SimInput.EStateData.StateTransitionCounts = zeros(SimInput.NumberStates,SimInput.NumberStates);
    
% If we are not using a prior equilibration run, then we need to run an
% equilibration simulation before we run the real simulation 
    if strcmp(SimInput.UsePriorEqRun,'n')
        
            SimInput.EStateData.InitEStatesGuess= EStates;
            
            % Run The Equilibration Simulation
            EQRun = 'y';
            disp('Starting Equilibration Run...')
            [SimOutput] = Run_Simulation_Engine(HexagonData,Options,SimInput,EQRun);
            
            % Plot the E state record
                if strcmp(Options.DisplayFigures,'y')
                    Show_EQ_And_Diagnostic_Results(Options,FigureHandles,HexagonData,SimOutput,EQRun);
                end
            
            % Set simulation output to be the input for the next
            % (non-equilibration) simulation run below 
                SimInput = SimOutput;
                SimInput.EStateData.EquilibratedEStates = SimInput.EStateData.EStatesRecord_Eq(:,:,end);
            
    % If we are using a prior equilibration run, then use that data to
    % start the real simulation 
    elseif strcmp(SimInput.UsePriorEqRun,'y')
        % Everything should already be loaded, so we just proceed to the post-equilibration run
        disp('Use Previous Equilibration Run')
    end    
   
    % Run The Real Simulation. Keep iterating until the minimum number of
    % virions fused has been reached, or until the maximum number of
    % iterations has been reached, as specified in the options
        EQRun = 'n';
        disp('Starting Simulation Run...')
        RunSimulationAgain = 'y';
        NumberIterations = 0;
        PreviousSimOutput = [];
        
        while strcmp(RunSimulationAgain,'y')
            NumberIterations = NumberIterations + 1;
            
            [SimOutput] = Run_Simulation_Engine(HexagonData,Options,SimInput,EQRun);

            % Aggregate the simulation out data across the different iterations
            [SimOutput] = Aggregate_Iteration_Data(PreviousSimOutput,SimOutput,NumberIterations,Options);
                        
            % Stop rerunning the simulation if the minimum number of fused
            % virions has been reached 
                if ~isempty(SimOutput.FusionWaitTimes)
                    NumberFusionEvents = length(SimOutput.FusionWaitTimes);
                else
                    NumberFusionEvents = 0;
                end
                disp(strcat('   Number virions fused =',num2str(NumberFusionEvents)));

                if NumberFusionEvents < Options.MinNumberFused
                    RunSimulationAgain = 'y';
                    disp('   Min number fused virions not reached, iterate simulation again...')
                else
                    RunSimulationAgain = 'n';
                end

            % Stop rerunning the simulation if the maximum number of iterations
            % has been reached, regardless of how many virions have fused 
                if NumberIterations >= Options.MaxSimIterations
                    disp('   Max number iterations reached, stop rerunning simulation.')
                    RunSimulationAgain = 'n';
                end
                
            % Record the simulation output for the next iteration
                PreviousSimOutput = SimOutput;
        end
    

    % Plot the E states
        if strcmp(Options.DisplayFigures,'y')
            Show_EQ_And_Diagnostic_Results(Options,FigureHandles,HexagonData,SimOutput,EQRun);
        end

end