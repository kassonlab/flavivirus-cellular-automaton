function [HexagonData] = Setup_Hexagon()

% ------------------------------
% This function generates a hexagon with a given number of monomers along 
% the top and along one of the sides. It automatically generates information 
% about which monomers pair together to form dimers, which trimers are 
% possible, and which trimers are considered nearby each other

% Output:
% HexagonData = data structure with the following fields:
    % MonomerInfoLibrary = Structure w/ information about each monomer, including 
        % the index of the dimer partner, the indexes of possible trimer 
        % partners, and the number of possible trimers in which the monomer could participate
    % DimerReferenceList = Vector to determine the index number of the dimer 
        % partner for a given monomer. For example, DimerReferenceList(4) 
        % would give you the index number of the dimer partner for monomer #4.
    % TrimerReferenceList = Matrix w/ the list of monomer index numbers in 
        % ascending order for each possible trimer. For example 
        % TrimerReferenceList(4) would give you the list of monomer index numbers 
        % of the monomers which participate in the trimer with trimer index number 4.
    % TrimerInfoLibrary = Structure w/ information about each possible trimer, 
        % ordered by trimer index number. This includes the index numbers of 
        % the monomers in the current trimer, list of nearby trimers, listed by 
        % their trimer index number, and the number of nearby possible trimers.
% ------------------------------
% A couple notes of explanation:

% Each monomer in the hexagon is given a set of coordinates (row number, 
% position number). These coordinates are organized a bit non-intuitively,
% so plot them out if it is strange to you. In retrospect, I would have set them up differently.

% Each monomer is also given an index number for easier reference 
% (instead of coordinates).

% Each trimer is given an index number for easier reference (so we can 
% refer to each unique trimer by its index number rather than a set of 
% 3 monomer index numbers).

% Written by Bob Rawle, 2017
% ------------------------------

% Initial options to configure:
SaveMatFile = 'y';
    NumberTop = 4;
    % The number of monomers along the top/bottom of the hexagon
    NumberSide = 4;
    % The number of monomers along the sides of the hexagon

% Define coordinates in the hexagon for each monomer
    MonomerIndexNumber = 0;

    for i = 1:NumberSide
        for j = 1:NumberTop + i - 1
            MonomerIndexNumber = MonomerIndexNumber +  1;
            MonomerInfoLibrary(MonomerIndexNumber).Coords = [i,j];
            MonomerCoords(MonomerIndexNumber,:) = [i,j];
        end
    end

    Cycle = 0;
    for i = NumberSide + 1: 2*NumberSide- 1
        Cycle = Cycle + 1;
        for j = 1:NumberTop + (2*NumberSide- 1) - i
            MonomerIndexNumber = MonomerIndexNumber +  1;
            MonomerInfoLibrary(MonomerIndexNumber).Coords = [i,j + Cycle];
            MonomerCoords(MonomerIndexNumber,:) = [i,j+ Cycle];
        end
    end

    NumberMonomers = MonomerIndexNumber;

% Define the dimers-note if you have an odd number of monomers, one of 
% them will end up not being attached to anybody (index number = 0)
    MonomerInfoLibrary(1).DimerIndexNumber = [];
    DimerReferenceList = zeros(1,NumberMonomers);

    for i = 1:NumberMonomers
        CurrentCoords = MonomerCoords(i,:);
        if isempty(MonomerInfoLibrary(i).DimerIndexNumber)
            PossibleBuddyCoords = CurrentCoords +[0,1];
            MonoIndexVector = MonomerCoords(:, 1) == PossibleBuddyCoords( 1, 1) & MonomerCoords(:, 2) == PossibleBuddyCoords( 1,2);
            PossibleBuddyIndex = find(MonoIndexVector);
            if length(PossibleBuddyIndex) ==1 && ...
                    isempty(MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber)

                MonomerInfoLibrary(i).DimerIndexNumber = PossibleBuddyIndex;
                MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber = i;
                DimerReferenceList(i) = PossibleBuddyIndex;
                DimerReferenceList(PossibleBuddyIndex) = i;
            else
                PossibleBuddyCoords = CurrentCoords +[1,1];
                MonoIndexVector = MonomerCoords(:, 1) == PossibleBuddyCoords( 1, 1) & MonomerCoords(:, 2) == PossibleBuddyCoords( 1,2);
                PossibleBuddyIndex = find(MonoIndexVector);
                if length(PossibleBuddyIndex) ==1 && ...
                    isempty(MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber)

                MonomerInfoLibrary(i).DimerIndexNumber = PossibleBuddyIndex;
                    MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber = i;
                    DimerReferenceList(i) = PossibleBuddyIndex;
                    DimerReferenceList(PossibleBuddyIndex) = i;
                else
                    PossibleBuddyCoords = CurrentCoords +[1,0];
                    MonoIndexVector = MonomerCoords(:, 1) == PossibleBuddyCoords( 1, 1) & MonomerCoords(:, 2) == PossibleBuddyCoords( 1,2);
                    PossibleBuddyIndex = find(MonoIndexVector);
                    if length(PossibleBuddyIndex) ==1 && ...
                        isempty(MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber)

                        MonomerInfoLibrary(i).DimerIndexNumber = PossibleBuddyIndex;
                        MonomerInfoLibrary(PossibleBuddyIndex).DimerIndexNumber = i;
                        DimerReferenceList(i) = PossibleBuddyIndex;
                        DimerReferenceList(PossibleBuddyIndex) = i;
                    else
                        MonomerInfoLibrary(i).DimerIndexNumber = 0;
                    end
                end
            end
        end
    end

% Define the possible trimers
    % TotalNumberofPossibleTrimersWDuplicates = 0;
    TrimerReferenceList = [];

    for i = 1:NumberMonomers
        CurrentNumberPossibleTrimers = 0;
        CurrentCoords = MonomerCoords(i,:);
        % Find all possible neighbors which could form a trimer (6 possibilities)
            PossibleN1 = CurrentCoords +[0,1];
            PossibleN2 = CurrentCoords +[1,1];
            PossibleN3 = CurrentCoords +[1,0];
            PossibleN4 = CurrentCoords +[0,-1];
            PossibleN5 = CurrentCoords +[-1,-1];
            PossibleN6 = CurrentCoords +[-1,0];

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN1( 1, 1) & MonomerCoords(:, 2) == PossibleN1(1,2);
            IndexN1 = find(MonoIndexVector);

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN2( 1, 1) & MonomerCoords(:, 2) == PossibleN2(1,2);
            IndexN2 = find(MonoIndexVector);

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN3( 1, 1) & MonomerCoords(:, 2) == PossibleN3(1,2);
            IndexN3 = find(MonoIndexVector);

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN4( 1, 1) & MonomerCoords(:, 2) == PossibleN4(1,2);
            IndexN4 = find(MonoIndexVector);

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN5( 1, 1) & MonomerCoords(:, 2) == PossibleN5(1,2);
            IndexN5 = find(MonoIndexVector);

            MonoIndexVector = MonomerCoords(:, 1) == PossibleN6( 1, 1) & MonomerCoords(:, 2) == PossibleN6(1,2);
            IndexN6 = find(MonoIndexVector);

        % Cycle through all possible trimers (6 possibilities if monomer is surrounded on all sides)
        if length(IndexN1) ==1 && length(IndexN2) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN1,IndexN2,CurrentNumberPossibleTrimers);
        end

        if length(IndexN2) ==1 && length(IndexN3) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN2,IndexN3,CurrentNumberPossibleTrimers);
        end

        if length(IndexN3) ==1 && length(IndexN4) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN3,IndexN4,CurrentNumberPossibleTrimers);
        end

        if length(IndexN4) ==1 && length(IndexN5) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN4,IndexN5,CurrentNumberPossibleTrimers);
        end

        if length(IndexN5) ==1 && length(IndexN6) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN5,IndexN6,CurrentNumberPossibleTrimers);
        end

        if length(IndexN6) ==1 && length(IndexN1) ==1
            CurrentNumberPossibleTrimers = CurrentNumberPossibleTrimers + 1;
            [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
                IndexN6,IndexN1,CurrentNumberPossibleTrimers);
        end

        MonomerInfoLibrary(i).NumberPossibleTrimers = CurrentNumberPossibleTrimers;

    % Diagnostics
    %     for j = 1:CurrentNumberPossibleTrimers
    %         TotalNumberofPossibleTrimersWDuplicates = TotalNumberofPossibleTrimersWDuplicates +1;
    %         TrimerIndicestoAdd = [i MonomerInfoLibrary(i).PossibleTrimerPairs(j,:)];
    %         TrimerList(TotalNumberofPossibleTrimersWDuplicates,:) = TrimerIndicestoAdd;
    %     end
    end

    %     TrimerReferenceList
    
% Determine which trimers (using the trimer index number) are considered nearby each other
    TrimerInfoLibrary = [];
    NumberUniqueTrimers = size(TrimerReferenceList, 1);

    for t = 1:NumberUniqueTrimers
        CurrTrimerMonoIndexes = TrimerReferenceList(t,:);
        CurrentTrimerIndex = t;
        TrimerInfoLibrary(CurrentTrimerIndex).MonoIndexes = CurrTrimerMonoIndexes;
        if CurrTrimerMonoIndexes(3) - 1 == CurrTrimerMonoIndexes(2)
            TrimerFaces = 'Down';
            DownIndex = CurrTrimerMonoIndexes(1);
                % Index of the monomer in the down position
        elseif CurrTrimerMonoIndexes(1) + 1 == CurrTrimerMonoIndexes(2)
            TrimerFaces = 'Up';
            UpIndex = CurrTrimerMonoIndexes(3);
                % Index of the monomer in the up position
        else
            ThrowError
        end

        % Depending on whether the trimer points up or down, search for all 6 
        % possibilities of nearby trimers and record the ones that exist in the
        % TrimerInfoLibrary. We search relative to the position of the monomer in the up or down position.
        if strcmp(TrimerFaces,'Up')

            % Define the relative coordinates of monomers in the possible nearby trimers
            NearbyTrimers(:,:,1) = [0,-1; 1, 0; 1,-1];
            NearbyTrimers(:,:,2) = [0, 1; 1, 1; 1, 2];
            NearbyTrimers(:,:,3) = [0, 1; 0, 2; -1, 1];
            NearbyTrimers(:,:,4) = [0, -1; 0, -2; -1, -2];
            NearbyTrimers(:,:,5) = [-2, -2; -2, -1; -3, -2];
            NearbyTrimers(:,:,6) = [-2, -1; -2, 0; -3, -1];

            UpCoords = MonomerCoords(UpIndex,:);
            NearbyTrimers(:,1,:) = NearbyTrimers(:,1,:) + UpCoords(1,1);
            NearbyTrimers(:,2,:) = NearbyTrimers(:,2,:) + UpCoords(1,2);

            [TrimerInfoLibrary] = Find_Nearby_Trimers(NearbyTrimers,MonomerCoords,TrimerReferenceList,...
                TrimerInfoLibrary,CurrentTrimerIndex);

        elseif strcmp(TrimerFaces,'Down')

            % Define the relative coordinates of monomers in the possible nearby trimers
            NearbyTrimers(:,:,1) = [0, 1; 0, 2; 1, 2];
            NearbyTrimers(:,:,2) = [0, 1; -1, 0; -1, 1];
            NearbyTrimers(:,:,3) = [0, -1; -1, -1; -1, -2];
            NearbyTrimers(:,:,4) = [0, -1; 0, -2; 1, -1];
            NearbyTrimers(:,:,5) = [2, 1; 2, 0; 3, 1];
            NearbyTrimers(:,:,6) = [2, 1; 2, 2; 3, 2];

            DownCoords = MonomerCoords(DownIndex,:);
            NearbyTrimers(:,1,:) = NearbyTrimers(:,1,:) + DownCoords(1,1);
            NearbyTrimers(:,2,:) = NearbyTrimers(:,2,:) + DownCoords(1,2);

            [TrimerInfoLibrary] = Find_Nearby_Trimers(NearbyTrimers,MonomerCoords,TrimerReferenceList,...
                TrimerInfoLibrary,CurrentTrimerIndex);

        end
    end
    
% Save the output data
HexagonData.MonomerInfoLibrary = MonomerInfoLibrary;
HexagonData.DimerReferenceList = DimerReferenceList;
HexagonData.TrimerReferenceList = TrimerReferenceList;
HexagonData.TrimerInfoLibrary = TrimerInfoLibrary;

if strcmp(SaveMatFile,'y')
    
    save(strcat('Hexagon_',num2str(NumberSide),'side_',num2str( NumberTop),'top.mat'),'HexagonData')
end
    
end

function [MonomerInfoLibrary,TrimerReferenceList] = Add_New_Trimer(i,MonomerInfoLibrary,TrimerReferenceList,...
    IndexN1,IndexN2,CurrentNumberPossibleTrimers)

    % Record the trimer pair for the current monomer
        MonomerInfoLibrary(i).PossibleTrimerPairs(CurrentNumberPossibleTrimers,:) = [IndexN1 IndexN2];
    
    % Search for the new trimer in the TrimerReferenceList
        NewTrimer = [i IndexN1 IndexN2];
        NewTrimer = sort(NewTrimer);
        
        if i > 1
            TrimerIndexVector = TrimerReferenceList(:, 1) == NewTrimer(1,1) &...
                TrimerReferenceList(:, 2) == NewTrimer(1,2) & TrimerReferenceList(:,3) == NewTrimer(1,3);
            TrimerIndex = find(TrimerIndexVector,1);

            % If the current trimer doesn't already exist in the reference 
            % list, add it to the TrimerReferenceList
            if isempty(TrimerIndex)
                TrimerReferenceList = [TrimerReferenceList; NewTrimer];
            end
        else
            TrimerReferenceList = [TrimerReferenceList; NewTrimer];
        end
end

function [TrimerInfoLibrary] = Find_Nearby_Trimers(NearbyTrimers,MonomerCoords,TrimerReferenceList,...
    TrimerInfoLibrary,CurrentTrimerIndex)

    NumberPossibilities = size(NearbyTrimers,3);
    NumberNearbyFound = 0;

    for n = 1:NumberPossibilities

        % First we convert the coordinates of the monomers in each 
        % of the nearby trimers to their respective monomer index number
            CurrNearbyTrimerMonoCoords = NearbyTrimers(:,:,n);
            for c =  1:3
                MonoIndexVector = MonomerCoords(:, 1) == CurrNearbyTrimerMonoCoords(c,1) & ...
                    MonomerCoords(:, 2) == CurrNearbyTrimerMonoCoords(c,2);
                MonoIndexFound = find(MonoIndexVector);
                if isempty(MonoIndexFound)
                    MonoIndexes = [];
                    break
                else 
                    MonoIndexes(c) = MonoIndexFound;
                end
            end

            MonoIndexes = sort(MonoIndexes);

        % If that nearby trimer exists, then we find the trimer index number of each set of monomers
        if length(MonoIndexes) == 3
            TrimerIndexVector = TrimerReferenceList(:, 1) == MonoIndexes(1,1) &...
                TrimerReferenceList(:, 2) == MonoIndexes(1,2) & ...
                TrimerReferenceList(:,3) == MonoIndexes(1,3);
            TrimerIndex = find(TrimerIndexVector);

            % If that trimer exists, record it as a nearby trimer
            if length(TrimerIndex) == 1
                NumberNearbyFound = NumberNearbyFound + 1;
                TrimerInfoLibrary(CurrentTrimerIndex).NearbyTrimers(NumberNearbyFound) = TrimerIndex;
                TrimerInfoLibrary(CurrentTrimerIndex).NumberNearbyTrimers = NumberNearbyFound;
            end
        end
    end
end