function [Options] = Setup_Options_Fit_CA

%-----Fit Options-----
Options.FitMethod = 'Constrained-PatternSearch';
    % Choose the fit method.
  
% Options.GuessArray = 'Large';
Options.GuessArray = [2000];
    % Your initial guess (these will be defined in Setup_Rate_Constants)
    % If you are scanning more than one guess, have additional 
    % guesses in separate rows. 'Large' means you will set up a larger guess
    % scan in Setup_Guess_Scan, at the bottom of this file.    
    if strcmp(Options.GuessArray,'Large')
        Options = Setup_Guess_Scan(Options);
    end
    
Options.LowerBounds =  [100];
Options.UpperBounds = [5000];
    % The lower and upper bounds on your fit parameters. Should be the same
    % size as your initial guess matrix. These values will only be used if
    % you are performing a constrained fit.

Options.LikelihoodMethod = 'KineticsAndEfficiencies';
% Options.LikelihoodMethod = 'Kinetics,GatedEfficiencies';
    % This determines which likelihood function will be used.
    % KineticsAndEfficiencies is what was used in the manuscript.

Options.KineticWeightConstant = 1;
    % Constant used in the log likelihood calculation to weight the
    % kinetic terms relative to the efficiencies. Set equal to 1 if you want
    % equal weighting. 

Options.NormalizeModelEfficiencies = 'n';
    % y or n. Normalize the efficiencies of the model to the first pH value
    % in the sequence when calculating the likelihood function. 

%-----Experimental data options-----
Options.ExpDataFilename = strcat(pwd,'/Data to fit/ZikaDataToFit.mat');
    % Path to experimental data file (.MAT) we will be fitting
Options.UseRelEfficiencies = 'y';
    % y or n. Whether or not to use relative experimental efficiencies,
    % normalizing to the value below
Options.EffMax = .191;
    % 0.191 = raw efficiency, Zika data pH 4.6

%------Simulation Options-----
Options.pHValues =  [4.6,5.5,5.8,6.1,6.6,6.9];
    % pH values you want to scan in the simulation
Options.TimeStep = 'Adapt';
    % Time per step, in seconds. Also, you can enter 'Adapt' instead of a number to
    % adaptively determine your time step, maximizing the time step to make
    % the fastest rate constant yield a probability of ProbCutoff, below.
Options.MinAdaptiveTimeStep =  0.001;
    % Minimum time step allowed if you are determining the time step adaptively.
Options.TotalTime = 340;
    % Total time of simulation, in seconds. Experimental data usually
    % collected for 340 seconds following pH drop. 
Options.NumberVirions = 1200;
    % Number of virions in each simulation
Options.MinNumberFused = 1000;
    % Minimum number of virions that need to fuse. If this number is not
    % reached within the first round of simulation using the number of virions
    % above, then multiple rounds will be performed and aggregated for the
    % given pH value until this minimum number is reached. Enter 0 if you
    % don't wish to have any minimum number.
Options.MaxSimIterations = 10;
    % Maximum number of simulation iterations allowed for a given pH value
    % to reach the minimum number of fused virions
Options.ProbCutoff = 0.1;
    % For a given transition, probabilities above this value will trigger a warning, results
    % may be unreliable. This value will also be used to determine
    % the time step if you select the option to determine it adaptively.

%-------Equilibration Options------
Options.NumberSteps_Eq = 100;
Options.TimeStep_Eq = 'Adapt';
    % Time step for the equilibration run, in seconds. Potentially, this
    % may be much longer than the time step for the real run, depending on
    % your rate constants. Also, you can enter 'Adapt' instead of a number to
    % adaptively determine your time step, maximizing the time step to make
    % the fastest rate constant yield a probability of ProbCutoff_Eq, below.
Options.EDataHowOften_Eq = floor(Options.NumberSteps_Eq/15);
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
    % don't want to use a preset file.
Options.DisplayFigures = 'y';
    % y or n. Whether or not to display any figures at all.
Options.CrudeParallel = 'y';
    % Choose y if you are running multiple instances of Matlab as a crude
    % parallelization. Make sure you set up your input properly into
    % Start_Fit_CA.
    
%----Legacy options, don't change----
Options.SaveEStateData = 'n';
    % y or n. Whether or not to save the state data for all E proteins in all viruses.
Options.EDataHowOften = 1;
    % How often to record the state data of all E proteins, in seconds.
Options.Diagnostics = 'n';

    %-----Gating test options/parameters---
    % See Calculate_Likelihood_Global_Fit.m for details on gate test.
    Options.GatePHValues = [4.6, 6.6];
%     Options.GatePHValues = [5, 6];
    % The 2 pH values whose efficiencies (from the model) will be compared to
    % see if they pass the gate test 
    Options.GateValue =  0.4;
    % The minimum difference in efficiencies needed to pass the gate test
    Options.HighEffValueCutoff = 0.8;
    % The minimum efficiency that the lowest pH value must have in order for the gate test to be valid

end

function Options = Setup_Guess_Scan(Options)
    
% Define individual parameters that you want to scan initial guesses
    Params{1,1} = [ 1, 500, 1300, 5000, 10000];
%     Params{1,2} = [ 1];
%     Params{1,3} = [ 0, 17, 25];

% Enumerate all possible combinations of your defined parameters. Note that
% the number of combinations can get quite large quite quickly...
    Options.GuessArray = allcomb(Params{1,1:size(Params,2)});    
end