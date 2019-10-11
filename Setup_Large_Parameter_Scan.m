function Options = Setup_Large_Parameter_Scan(Options,SaveFolderDir)
    
% Define individual parameters that you want to scan
    Params{1,1} = [100000,250000]; %[40000];
    Params{1,2} = [6.8]; 
    Params{1,3} = [5.6];
    Params{1,4} = [100,25000,50000];
    Params{1,5} = [15000,100000,500000];
    Params{1,6} = [.05,.15,.4];

% Enumerate all possible combinations of your defined parameters. Note that
% the number of combinations can get quite large quite quickly...
    Options.ParameterArray = allcomb(Params{1,1:size(Params,2)});

    % Save parameter reference file
    if strcmp(Options.SaveData,'y')
        ParameterReference.ParameterArray = Options.ParameterArray;
        ParameterReference.Params = Params;
        save(strcat(SaveFolderDir,'/',Options.FileLabel,';ParameterReference','.mat'),'ParameterReference')
    end
    
end