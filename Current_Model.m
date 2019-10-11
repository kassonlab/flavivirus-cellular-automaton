function [SimInput] = Current_Model(CurrentParameters,SimInput)
    % Hardcoded states for each E protein monomer:
        %1: retracted (starting state)
        %2: extended
        %3: in a trimer
        %4: hemifused
    % You can add other states and create connections between them as you
    % see fit, but don't modify the states above as they are hardcoded into
    % the simulation  
    
    % Define number of states in the model
    SimInput.PossibleStates = 1:4;
    SimInput.NumberStates = length(SimInput.PossibleStates);

    % Extract pH value
    pH = SimInput.pH;
        HConc = 10^-pH;
    
    % Set up rate constant matrix for the target pH. k(starting state, ending state).
    % All rate constants have time units of seconds^-1. Those multiplied by
    % proton concentration also have units of M^-1.
    % Note: this matrix is the transposition of how a system of
    % differential equations is usually set up, but I think that setting it
    % up this way is more intuitive, which is why it is defined as above.
        k = zeros(SimInput.NumberStates,SimInput.NumberStates);
        
        k23_base = CurrentParameters(1); %5000
        k12_base = 100*k23_base;
            pK_12 = 6.8;
            pK_23 = 8;

        k(1,2) = k12_base * HConc;
        k(2,1)= k12_base * 10^-pK_12;
        k(2,3) = k23_base * HConc;
        k(3,2) = 1*10^-pK_23;
        k(3,4) = 3e-3;
        
    % Set up rate constant matrix for equilibration run (pH 7.4).
        k_Eq = k;
            pH_Eq = 7.4;
%         k_Eq(1,2) = k12_base * 10^-pH_Eq;
        k_Eq(1,2) = 0;

        k_Eq(2,3) = 0;
        k_Eq(3,2) = 0;
        k_Eq(3,4) = 0;
        
    % Record rate constant matrices as simulation inputs
    SimInput.RateConstantMatrix = k;
    SimInput.RateConstantMatrix_Eq = k_Eq;
    
	% Define other parameters which will be used in the simulation
    SimInput.CoopFactor = 1;
        % Factor to increase the rate constant of a given monomer becoming extended
        % if its dimer pair is already out
    SimInput.MinNumTrimersForFusion = 2;
        % Minimum number of adjacent trimers required for fusion.

end