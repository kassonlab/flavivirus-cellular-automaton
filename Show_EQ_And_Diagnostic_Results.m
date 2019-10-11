function Show_EQ_And_Diagnostic_Results(Options,FigureHandles,HexagonData,SimOutput,EQRun)
    
% Display the results of the equilibration run to make sure that it has actually equilibrated
if strcmp(EQRun,'y')
    set(0,'CurrentFigure',FigureHandles.EQWindow)
    clf

    NumberMonomers = length(HexagonData.MonomerInfoLibrary);
    EData = SimOutput.EStateData.EStatesRecord_Eq;
    ETime = SimOutput.EStateData.TimeVector_Eq;
    
    for State = SimOutput.PossibleStates
        StateData = floor(EData) == State;
        StateSumData(State,1:length(ETime)) = sum(sum(StateData))./(Options.NumberVirions*NumberMonomers);

        set(0,'CurrentFigure',FigureHandles.EQWindow)
        hold on
        plot(ETime,StateSumData(State,:))
        drawnow

        DLegend{1,State} = strcat('EState=',num2str(State));
    end

    set(0,'CurrentFigure',FigureHandles.EQWindow)
    legend(DLegend,'Location','best');
    xlabel('Time');
    ylabel('EState Frequency');
    title('Equilibration Run')

    % Display the diagnostic results, if selected
    if strcmp(Options.Diagnostics,'y')
        % Show the transition counts
        disp('Equilibration Results:')
        disp('   E Transition (Equilibration) Counts = ')
        disp(SimOutput.EStateData.StateTransitionCounts_Eq);
        disp('      - - - - - - - - - - - - - - ')
    end
else 
    % Display the diagnostic results, if selected
    if strcmp(Options.Diagnostics,'y')
        set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        clf

        NumberMonomers = length(HexagonData.MonomerInfoLibrary);
        EData = SimOutput.EStateData.EStatesRecord;
        ETime = SimOutput.EStateData.TimeVector;
        for State = SimOutput.PossibleStates
            StateData = floor(EData) == State;
            StateSumData(State,1:length(ETime)) = sum(sum(StateData))./(Options.NumberVirions*NumberMonomers);

            set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
            hold on
            plot(ETime,StateSumData(State,:))
            drawnow

            DLegend{1,State} = strcat('EState=',num2str(State));
        end

        set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        legend(DLegend,'Location','best');
        xlabel('Time');
        ylabel('EState Frequency');
        
        % Show the transition counts
            disp('Simulation Results:')
            disp('E Transition Counts = ')
            disp(SimOutput.EStateData.StateTransitionCounts);
            disp('---------------------------------')
    end
end

end