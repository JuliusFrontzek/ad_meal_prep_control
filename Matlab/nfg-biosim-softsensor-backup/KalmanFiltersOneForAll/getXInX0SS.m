%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  [xIn,x0SS] = getXInX0SS(flagModel,flagFrac,resultsSoeren,varargin)
% returns the inlet concentrations xIn and the initial concentrations x0SS
% for computing the state trajectory into steady state 

% xIn - inlet concentrations [g/l]
% x0SS - initial condition for state trajectory into steady state (from ...
% where in turn the dynamic simulation starts) 
% flagModel -       3: ADM1-R3; 4: ADM1-R4
% flagFrac -        0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% resultsSoeren -   struct from Sörens GitHub containing both initial ...
% conditions and inlet concentrations
% optional: fracChFast - share of fast carbohydrates in total carbohydrates

if flagFrac == 1
    fracChFast = varargin{1}; 
end

%% inlet concentrations XIn:
% miscellaneous for inlet concentrations:
xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
nStatesSoeren = length(xInPre); % Soeren has no Ash and no 2nd CH-fraction
xAshIn = 14; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

switch flagModel
    case 3
        % beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
        % fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
        % der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
        % ist es folgendermaßen:
        ScatIn = xInPre(11);
        SanIn = xInPre(12); 
        SionIN = ScatIn - SanIn; 
        % adapt inlet concentrations for slightly different state indexing in
        % Simon's ADM1-R3 models (X_ash, S_ion = S_cat - S_an):
        if flagFrac == 0
            xIn = nan(nStatesSoeren,1); % allocate memory
            xIn([1:10,13:end]) = xInPre([1:10,13:end])'; % all states except Xash and Sion
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            xIn(11) = xAshIn;  
            xIn(12) = SionIN; 
        else
            % consider second CH fraction:
            xIn = nan(nStatesSoeren + 1,1); % allocate memory
            xIn([1:6,8:11,14:end]) = xInPre([1:10,13:end])'; % all states except Xchslow, Xash and Sion
            xIn(7) = 0;     % all inlet carbohydrates are assigned to the 
            % fast fraction. Separation in fast/slow only in ODE function
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            xIn(12) = xAshIn;  
            xIn(13) = SionIN; 
        end
    case 4
        if flagFrac == 0
            xIn = nan(nStatesSoeren + 1,1); % mind that Soeren has no Ash!   
            xIn([1:8,10:end]) = xInPre'; 
            xIn(9) = xAshIn;
        else 
            xIn = nan(nStatesSoeren + 2,1); % mind that Soeren has no Ash and no 2nd CH fraction   
            xIn([1:5,7:9,11:end]) = xInPre'; 
            xIn(6) = 0;     % assign all carbohydrates in to CH_fast. 
            % only separate total carbs in in convection terms of ODE!
            xIn(10) = xAshIn;
        end
end

%% initial condition x0SS for transition into steady state:
% miscellaneous for initial condition: 
x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
xAsh0 = 1;                  % selbst gewählt

switch flagModel
    case 3
        Scat0 = x0Pre(11); 
        San0 = x0Pre(12); 
        Sion0 = Scat0 - San0; 
        % adapt initial condition for slightly different state indexing in standard 
        % control notation (X_ash, S_ion = S_cat - S_an): 
        if flagFrac == 0
            x0SS = nan(nStatesSoeren,1); 
            x0SS(1:10) = x0Pre(1:10);       % all states up to X_ac 
            x0SS(13:end) = x0Pre(13:end);   % all other states 
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            x0SS(11) = xAsh0;
            x0SS(12) = Sion0;   
            % since the dissociation constants KaIN and Kaco2 were increased 
            % compared with Sörens false values on GitHub, increase the
            % initial values of biomass to sustain a positive steady state:
            x0SS(9) = 1.5*x0SS(9);      % microbial biomass 
            x0SS(10) = 1.5*x0SS(10);    % methanogens
        else
            % consider additional CH fraction
            % XY: Werte für XchS hinzufügen!
            x0SS = nan(nStatesSoeren + 1,1); % allocate memory
            x0SS([1:5,8:11,14:end]) = xInPre([1:5,7:10,13:end])'; % all states except Xch (slow and fast), Xash and Sion
            % overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
            x0SS(12) = xAsh0;
            x0SS(13) = Sion0;
            % add corbohydrate inlet concentrations: 
            Xch0 = x0Pre(6);   
            x0SS(6) = fracChFast*Xch0; 
            x0SS(7) = (1-fracChFast)*Xch0;
        end

    case 4
        if flagFrac == 0
            x0SS = nan(nStatesSoeren + 1,1);    % mind that Soeren has no ash! 
            x0SS([1:8,10:end]) = x0Pre;         % all states up to X_ash 
            x0SS(9) = xAsh0; 
        else 
            % XY: Werte für XchS einfügen!
            x0SS = nan(nStatesSoeren + 2,1); % mind that Soeren has no Ash and no 2nd CH fraction   
%             x0SS = [0.091, 0.508, 0.944, 956.97, fracChFast*x0ch, (1-fracChFast)*x0ch, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 
            x0chTot = x0Pre(5); 
            x0SS([1:4,7:9,11:end]) = x0Pre([1:4,6:end]); 
            x0chFast = fracChFast*x0chTot; 
            x0chSlow = (1 - fracChFast)*x0chTot; 
            x0SS(5:6) = [x0chFast;x0chSlow]; 
            x0SS(10) = xAsh0; 
        end

end

end