function [SimOutput] = Aggregate_Iteration_Data(OldData,NewData,NumberIterations,Options)

SimOutput = NewData;

    if ~isempty(OldData)
       SimOutput.FusionWaitTimes = [OldData.FusionWaitTimes,NewData.FusionWaitTimes];
       SimOutput.EStateData.EStatesRecord = cat(1,OldData.EStateData.EStatesRecord,NewData.EStateData.EStatesRecord);
       SimOutput.EStateData.StateTransitionCounts = ...
           OldData.EStateData.StateTransitionCounts + NewData.EStateData.StateTransitionCounts;
    end
    
SimOutput.NumberVirions = NumberIterations*Options.NumberVirions;
end