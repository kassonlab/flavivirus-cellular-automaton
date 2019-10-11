function [Options] = Setup_Options()

Options.SaveData = 'y';
    % y or n, whether or not to save the simulation data as a .mat file
Options.FileLabel = 'test';
    % File label for the saved .mat file (be careful you don't overwrite anything)

%------Simulation Options-----
Options.pHValues = [5.0,5.25,5.5,5.75,6.0,6.25];
    % pH values you want to scan in the simulation
Options.ParameterArray = 'Large';
    % Your scan parameter (this will be pH or a specific rate constant 
    % or something). If you are scanning more than one parameter, have additional 
    % parameters in separate rows. NaN indicates you are not varying
    % anything. 'Large' means you will set up a larger parameter scan in Setup_Large_Parameter_Scan
Options.TimeStep = 'Adapt';
    % Time per step, in seconds. Also, you can enter 'Adapt' instead of a number to
    % adaptively determine your time step, maximizing the time step to make
    % the fastest rate constant yield a probability of ProbCutoff, below.
Options.MinAdaptiveTimeStep =  .001;
    % Minimum time step allowed if you are determining the time step adaptively.
Options.TotalTime = 340;
    % Total time of simulation, in seconds. Experimental data usually
    % collected for 340 seconds following pH drop. 
Options.NumberVirions = 500;
    % Number of virions in each simulation
Options.MinNumberFused = 200;
    % Minimum number of virions that need to fuse. If this number is not
    % reached within the first round of simulation using the number of virions
    % above, then multiple rounds will be performed and aggregated for the
    % given pH value until this minimum number is reached. Enter 0 if you
    % don't wish to have any minimum number.
Options.MaxSimIterations = 3;
    % Maximum number of simulation iterations allowed for a given pH value
    % to reach the minimum number of fused virions
Options.SaveEStateData = 'y';
    % y or n. Whether or not to save the state data for all E proteins in all viruses.
Options.EDataHowOften = 1;
    % How often to record the state data of all E proteins, in seconds.
Options.ProbCutoff = 0.1;
    % For a given transition, probabilities above this value will trigger a warning, results
    % may be unreliable. This value will also be used to determine
    % the time step if you select the option to determine it adaptively.
    
%-------Equilibration Options------
Options.NumberSteps_Eq = 3000;
Options.TimeStep_Eq = 0.1; %'Adapt'
    % Time step for the equilibration run, in seconds. Potentially, this
    % may be much longer than the time step for the real run, depending on
    % your rate constants. Also, you can enter 'Adapt' instead of a number to
    % adaptively determine your time step, maximizing the time step to make
    % the fastest rate constant yield a probability of ProbCutoff, below.
Options.EDataHowOften_Eq = floor(Options.NumberSteps_Eq/10);
    % How often to record the state data of all E proteins, in number of steps.
Options.ReuseEq = 'y';
    % y or n. Whether or not to re-use the endpoint of the equilibration run of
    % the first pH value within a given parameter set as the equilibrated state (e.g. starting state) of
    % all subsequent pH values. If you choose n, a separate equilibration run
    % will be performed for each pH value, increasing computational time.   
Options.ProbCutoff_Eq = 0.1;
    % Probabilities above this value will trigger a warning, results
    % may be unreliable. This value will also be used to determine
    % the time step if you select the option to determine it adaptively to
    % your rate constants.

%----Other Options----
Options.PresetHexagonFilename = strcat(pwd,'/Hexagon Choices/Hexagon_4side_3top.mat');
% Options.PresetHexagonFilename = [];
    % The name of the .mat file name (should be in the working directory) of 
    % a hexagon data file which has already been compiled. Leave empty if you 
    % don't want to use a preset file. In this case, a hexagon will be set
    % up according to Setup_Hexagon.m. Make sure that file is properly
    % configured.  
Options.DisplayFigures = 'n';
    % y or n. Whether or not to display any figures at all.
Options.Diagnostics = 'n';
    % This will show diagnostic plot(s), including how the E states vary
    % with time, as well as some diagnostic metrics printed out to the
    % command line
Options.Parallel = 'n';
    % Choose y if you are running multiple instances of Matlab as a crude
    % parallelization. This will split up the parameter scan according to your
    % input into Start_Zika_Cellular_Automaton.  
end