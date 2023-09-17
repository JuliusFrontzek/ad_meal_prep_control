% Funktion zur Multimodell Trajektorienplanung

function optTP = MMTP(TP,w_aic)

u0 = TP.u(2,:);

%% Aufgabe 2b:
optTP = fmincon(@guete_TP,u0,[],[],[],[],TP.CONS.umin,TP.CONS.umax,[],[]);

%% Gütefunktion für MMTP

    function I = guete_TP(uTP)

        TP.u = [TP.u(1,:); uTP; TP.u(3,:)];       
        p = TP.p;
        for idxMess = 1:length(TP)

        %% Simulation

        x0 = TP(idxMess).x0;
        x0_ = x0;

        for Modell = 1:3

            for idxT = 2:length(TP(idxMess).t)
            tspan = [TP(idxMess).t(idxT-1) TP(idxMess).t(idxT)];

                if Modell == 1
                    [T,X1] = ode45(@biomodell_1_zdgl,tspan,x0,[],TP(idxMess).u,p(Modell).p);
                    x0 = X1(end,:)';
                elseif Modell == 2
                    [T,X2] = ode45(@biomodell_2_zdgl,tspan,x0,[],TP(idxMess).u,p(Modell).p);
                    x0 = X2(end,:)';
                elseif Modell == 3
                    [T,X3] = ode45(@biomodell_3_zdgl,tspan,x0,[],TP(idxMess).u,p(Modell).p);
                    x0 = X3(end,:)';
                end

            end

        end % idxT
        X_end = [X1(end,1) X2(end,1) X3(end,1)];

        %% 2b: Gütefunktional: gewichtete Summe der letzten Werte der Biomasse
        I = 

        end % for idxMess		

    end

end
