%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  [xIn,x0SS] = getXInX0SS(flagModel,flagFrac,resultsSoeren)
% XY: Zweck der Funktion und Argumente beschreiben

% inlet concentrations XIn:
switch flagModel
    case 3
        if flagFrac == 0
            xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
            % beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
            % fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
            % der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
            % ist es folgendermaßen:
            ScatIn = xInPre(11);
            SanIn = xInPre(12); 
            SionIN = ScatIn - SanIn; 
            % adapt inlet concentrations for slightly different state indexing in
            % Simon's ADM1-R3 model (X_ash, S_ion = S_cat - S_an): 
            nStates = length(xInPre); 
            xIn = nan(nStates,1); % allocate memory
            xIn([1:10,13:end]) = xInPre([1:10,13:end])'; % all states except Xash and Sion
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            xAshIn = 14; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)
            xIn(11) = xAshIn;  
            xIn(12) = SionIN; 
        else
            % Gleichungen für R3-frac nachtragen
        end
end

% initial condition x0SS for transition into steady state:
switch flagModel
    case 3
        x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
        Scat0 = x0Pre(11); 
        San0 = x0Pre(12); 
        Sion0 = Scat0 - San0; 
        % adapt initial condition for slightly different state indexing in standard 
        % control notation (X_ash, S_ion = S_cat - S_an): 
        x0SS = nan(nStates,1); 
        x0SS(1:10) = x0Pre(1:10);       % all states up to X_ac 
        x0SS(13:end) = x0Pre(13:end);   % all other states 
        % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
        xAsh0 = 1; % assume 1 g/l initial ash concentration
        x0SS(11) = xAsh0;
        x0SS(12) = Sion0;   
        % since the dissociation constants KaIN and Kaco2 were increased, increase
        % the initial values of biomass to sustain a positive steady state:
        x0SS(9) = 1.5*x0SS(9);      % microbial biomass 
        x0SS(10) = 1.5*x0SS(10);    % methanogens
    
end

end