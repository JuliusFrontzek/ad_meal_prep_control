	function [FM,EW,invFM] = guete_ovp_FIM(VP,uVP,FM_old,invC)
		VP.u = [VP.u(1,:); uVP; VP.u(3,:)];
		
%% -- Fishermatrix
		
		FM = biomodell_FIM(VP,VP.p,invC) + FM_old;
		
%% -- A-Kriterium
		B=eye(length(VP.p));
		%invFM = inv(FM);
        invFM=FM\B;
		EW = eig(invFM);
		I = sum(EW);