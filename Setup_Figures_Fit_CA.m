function [FigureHandles] = Setup_Figures_Fit_CA(Options)

%     if strcmp(Options.Diagnostics,'y')
%             FigureHandles.DiagnosticWindow = figure(3);
%             set(FigureHandles.DiagnosticWindow,'Position',[848   364   450   341]);
%             clf
%     end
    
%         FigureHandles.BackgroundTraceWindow = figure(4);
%         set(FigureHandles.BackgroundTraceWindow,'Position',[452 -130 450 341]);
if strcmp(Options.DisplayFigures,'y')
        Width = 400;
        Height = 300;
        FigureHandles.ExpPlot = figure(1);
        set(FigureHandles.ExpPlot,'Position',[6   479   Width  Height]);
        clf
        
        FigureHandles.ExpNormalizedPlot= figure(2);
        set(FigureHandles.ExpNormalizedPlot,'Position',[407 405   Width  Height]);
        clf
        
        FigureHandles.ModelPlot = figure(3);
        set(FigureHandles.ModelPlot,'Position',[6  50   Width  Height]);
        clf
        
        FigureHandles.ModelNormalizedPlot = figure(4);
        set(FigureHandles.ModelNormalizedPlot,'Position',[407 55   Width  Height]);
        clf
        
        FigureHandles.EfficiencyPlot = figure(5);
        set(FigureHandles.EfficiencyPlot,'Position',[808 405   Width  Height]);
        clf
else
    FigureHandles = [];
end
end