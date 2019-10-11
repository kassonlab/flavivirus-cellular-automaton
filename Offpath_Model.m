function [SimInput] = Offpath_Model(CurrentParameters,SimInput)
    % Hardcoded states for each E protein monomer:
        %1: retracted (starting state)
        %2: extended
        %3: in a trimer
        %4: hemifused
    % You can add other states and create connections between them as you
    % see fit, but don't modify the states above as they are hardcoded into
    % the simulation  
    
    % Additional states added:
        %5 off-pathway state
        
    SimInput.PossibleStates = 1:5;
    SimInput.NumberStates = length(SimInput.PossibleStates);

    pH = SimInput.pH;
        HConc = 10^-pH;
    
    % Set up rate constant matrix for the target pH. k(starting state, ending state).
    % All rate constants have time units of seconds^-1. Those multiplied by
    % proton concentration also have units of M^-1.
        k = zeros(SimInput.NumberStates,SimInput.NumberStates);
        
        k23_base = CurrentParameters(1);
        k12_base = CurrentParameters(5); %k23_base*100;
        k25_base = k23_base * CurrentParameters(6);
            pK_12 = CurrentParameters(2);
            pK_25 = CurrentParameters(3);
           %pK_23 = CurrentParameters(4);
            
        k(2,5) = k25_base *10^-pK_25;
%        k(5,2) = k25_base * HConc;
        
        k(1,2) = k12_base * HConc; 
        k(2,1)= k12_base * 10^-pK_12;
        
        k(2,3) = k23_base * HConc;
%        k(3,2) = k23_base * 10^-pK_23;
        
        k(3,4) = .0047;
        
    % Set up rate constant matrix for equilibration run (pH 7.4).
        k_Eq = k;
            pH_Eq = 7.4;
            
%        k_Eq(2,5) = 0;
        %k_Eq(5,2) = 0;
        
        k_Eq(1,2) = k12_base * 10^-pH_Eq;
        %k_Eq(1,2)=0;
        %k_Eq(2,1)=0;
        
       k_Eq(2,3) = k23_base * 10^-pH_Eq;
%         k_Eq(2,3) = 0;
%         k_Eq(3,2) = 0;
%         
%         k_Eq(3,4) = 0;
        
    % Record rate constant matrices as simulation inputs
    SimInput.RateConstantMatrix = k;
    SimInput.RateConstantMatrix_Eq = k_Eq;
    
	% Define other parameters which will be used in the simulation
    SimInput.CoopFactor = CurrentParameters(4);
        % Factor by which probability of a given monomer becoming active 
        % will increase if its dimer pair is already out
    SimInput.MinNumTrimersForFusion = 2;
        % Minimum number of adjacent trimers required for fusion.

end
    