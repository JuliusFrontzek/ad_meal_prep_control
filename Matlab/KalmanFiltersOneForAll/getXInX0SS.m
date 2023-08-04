%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  [xIn,x0SS] = getXInX0SS(flagModel,flagFrac,resultsSoeren)
% returns the inlet concentrations xIn and the initial concentrations x0SS
% for computing the state trajectory into steady state 

% xIn - inlet concentrations [g/l]
% x0SS - initial condition for state trajectory into steady state (from ...
% where in turn the dynamic simulation starts) 
% flagModel -       3: ADM1-R3; 4: ADM1-R4
% flagFrac -        0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% resultsSoeren -   struct from Sörens GitHub containing both initial ...
% conditions and inlet concentrations

%% inlet concentrations XIn:
% miscellaneous for inlet concentrations:
xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
nStatesSoeren = length(xInPre); % Soeren has no Ash and no 2nd CH-fraction
xAshIn = 14; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

switch flagModel
    case 3
        if flagFrac == 0
            % beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
            % fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
            % der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
            % ist es folgendermaßen:
            ScatIn = xInPre(11);
            SanIn = xInPre(12); 
            SionIN = ScatIn - SanIn; 
            % adapt inlet concentrations for slightly different state indexing in
            % Simon's ADM1-R3 model (X_ash, S_ion = S_cat - S_an):
            xIn = nan(nStatesSoeren,1); % allocate memory
            xIn([1:10,13:end]) = xInPre([1:10,13:end])'; % all states except Xash and Sion
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            xIn(11) = xAshIn;  
            xIn(12) = SionIN; 
        else
            % Gleichungen für R3-frac nachtragen
            % XY: du wirst dann fracChFast als Argument in der Funktion
            % getXInX0SS.m brauchen!
        end
    case 4
        if flagFrac == 0
            %      S_ch4, S_IC,S_IN,  S_h2o,   X_ch,   X_pr,  X_li,  X_bac,X_ash,  S_ch4,g, S_co2,g
    %       xIn = [0,     0,   0.592, 960.512, 23.398, 4.75,  1.381, 0,    17,     0,    0]'; % [g/l], 
            xIn = nan(nStatesSoeren + 1,1); % mind that Soeren has no Ash!   
            xIn([1:8,10:end]) = xInPre'; 
            xIn(9) = xAshIn;
        else 
            % XY: Eingangs-Konzentrationen für XchS einfügen mit
            % fracChFast!
        end
end

%% initial condition x0SS for transition into steady state:
% miscellaneous for initial condition: 
x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
xAsh0 = 1;                  % selbst gewählt

switch flagModel
    case 3
        if flagFrac == 0
            Scat0 = x0Pre(11); 
            San0 = x0Pre(12); 
            Sion0 = Scat0 - San0; 
            % adapt initial condition for slightly different state indexing in standard 
            % control notation (X_ash, S_ion = S_cat - S_an): 
            x0SS = nan(nStatesSoeren,1); 
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
        else
            % XY: Werte für XchS hinzufügen!
        
        end

    case 4
        if flagFrac == 0
            x0SS = nan(nStatesSoeren + 1,1);    % mind that Soeren has no ash! 
            x0SS([1:8,10:end]) = x0Pre;         % all states up to X_ash 
            x0SS(9) = xAsh0; 
        else 
            % XY: Werte für XchS einfügen!

        end

end

end