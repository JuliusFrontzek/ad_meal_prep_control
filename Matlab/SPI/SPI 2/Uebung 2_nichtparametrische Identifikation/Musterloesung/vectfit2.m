function [SER,poles,rmserr,fit]=vectfit2(f,s,poles,weight,VF);
%
%function [SER,poles,rmserr,fit]=vectfit2(f,s,poles,weight,VF)
% 
%     ===========================================================
%     =   Vector Fitting                                        =
%     =   Version 2.1                                           =
%     =   Last revised: 27.10.2005                              = 
%     =   Bjorn Gustavsen                                       =
%     =   SINTEF Energy Research, N-7465 Trondheim, NORWAY      =
%     =   bjorn.gustavsen@sintef.no                             =
%     =   http://www.energy.sintef.no/Produkt/VECTFIT/index.asp =
%     ===========================================================
%
% PURPOSE : Approximate f(s) with a state-space model 
%
%         f(s)=C*(s*I-A)^(-1)*B +D +s*E
%                       
%           where f(s) is a singe element or a vector of elements.  
%           When f(s) is a vector, all elements become fitted with a common
%           pole set.
%
% INPUT :
%
% f(s) : function (vector) to be fitted. 
%        dimension : (Nc,Ns)  
%                     Nc : number of elements in vector
%                     Ns : number of frequency samples 
% 
% s : vector of frequency points [rad/sec] 
%        dimension : (1,Ns)  
%
% poles : vector of initial poles [rad/sec]
%         dimension : (1,N)  
%
% weight: the rows in the system matrix are weighted using this array. Can be used 
%         for achieving higher accuracy at desired frequency samples. 
%         If no weighting is desired, use: weight=ones(1,Ns). 
%
%         Two dimensions are allowed:
%           dimension : (1,Ns) --> Common weighting for all vector elements.   
%           dimension : (Nc,Ns)--> Individual weighting for vector elements.  
%
% VF.relax==1 --> Use relaxed nontriviality constraint 
% VF.relax==0 --> Use nontriviality constraint of "standard" vector fitting
%
% VF.kill=0 --> unstable poles are kept unchanged
% VF.kill=1 --> unstable poles are deleted
% VF.kill=2 --> unstable poles are 'flipped' into the left half plane
%            (kill=2 is the recommended choice)
%
% VF.asymp=1 --> Fitting with D=0,  E=0 
% VF.asymp=2 --> Fitting with D~=0, E=0 
% VF.asymp=3 --> Fitting with D~=0, E~=0 
%
% VF.spy1=1 --> Plotting, after pole identification (A)
%              figure(3): magnitude functions
%                cyan trace  : (sigma*f)fit              
%                red trace   : (sigma)fit
%                green trace : f*(sigma)fit - (sigma*f)fit
%
% VF.spy2=1 --> Plotting, after residue identification (C,D,E) 
%              figure(1): magnitude functions
%              figure(2): phase angles  
%
% VF.logx=1 --> Plotting using logarithmic absissa axis             
%
% VF.logy=1 --> Plotting using logarithmic ordinate axis
%
% VF.errplot=1   --> Include deviation in magnitude plot
%
% VF.phaseplot=1 -->Show plot also for phase angle
%
%
% VF.skip_pole=1 --> The pole identification part is skipped, i.e (C,D,E) 
%                    are identified using the initial poles (A) as final poles.
%                 
% VF.skip_res =1 --> The residue identification part is skipped, i.e. only the 
%                    poles (A) are identified while C,D,E are returned as zero. 
%
% VF.use_normal=1 -->Solving Least Squares (LS) systems using the Normal Equations.  
%              =0 -->Solving LS systems using QR decomposition 
%
% VF.use_sparse=1 -->Sparsity is used for formulating and solving LS system of the 
%                    pole identification problem. 
%              =0 -->Full arithmetic is used (not very useful).  
%
% VF.cmplx_ss  =1 -->The returned state-space model has real and complex conjugate 
%                    parameters. Output variable A is diagonal (and sparse). 
%              =0 -->The returned state-space model has real parameters only.
%                    Output variable A is square with 2x2 blocks (and sparse).
%
% OUTPUT :
% 
%     fit(s) = C*(s*I-(A)^(-1)*B +D +s.*E
%
% SER.A(N,N)    : A-matrix (sparse). If cmplx_ss==1: Diagonal and complex. 
%                                Otherwise, square and real with 2x2 blocks. 
%                           
% SER.B(N,1)    : B-matrix. If cmplx_ss=1: Column of 1's. 
%                       If cmplx_ss=0: contains 0's, 1's and 2's)
% SER.C(Nc,N)   : C-matrix. If cmplx_ss=1: complex
%                       If cmplx_ss=0: real-only
% SERD.D(Nc,1)  : constant term (real). Is non-zero if asymp=2 or 3.
% SERE.E(Nc,1)  : proportional term (real). Is non-zero if asymp=3.
%
% poles(1,N)    : new poles 
%
% rmserr(1) : root-mean-square error of approximation for f(s). 
%                   (0 is returned if skip_res==1)  
% fit(Nc,Ns): Rational approximation at samples. (0 is returned if
%             skip_res==1).
%
%
% APPROACH: The identification is calculated using Vector Fitting, see
% reference in text box below. A modification has been introduced in v2.0 which 
% makes convergence substantially faster and more reliable than in v1.0, see
% B. Gustavsen, "Improving the pole relocating properties of Vector
% Fitting", IEEE Trans. Power Delivery, accepted.
%
%
%******************************************************************************** 
% NOTE: This program is in the public domain and may be used by anyone. If the  *
%       program code (or a modified version) is used in a scientific work, or   *
%       in a commercial program, then reference should be made to the following:*
%       B. Gustavsen and A. Semlyen, "Rational approximation of frequency       *
%       domain responses by Vector Fitting", IEEE Trans. Power Delivery,        *
%       vol. 14, no. 3, pp. 1052-1061, July 1999.                               * 
%********************************************************************************

%Tolerances used by relaxed version of vector fitting
TOLlow=1e-8; TOLhigh=1e8;

[a,b]=size(poles); 
if s(1)==0 && a==1
  if poles(1)==0 && poles(2)~=0 
    poles(1)=-1;
  elseif poles(2)==0 && poles(1)~=0 
    poles(2)=-1;
  elseif poles(1)==0 && poles(2)==0 
    poles(1)=-1+i*10; poles(2)=-1-i*10;
  end
end  

if (VF.relax~=0) && (VF.relax)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.relax:  ' num2str(VF.asymp)]),return
end
if (VF.asymp~=1) && (VF.asymp)~=2 && (VF.asymp)~=3
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.asymp:  ' num2str(VF.asymp)]),return
end
if (VF.kill~=0) && (VF.kill~=1) && (VF.kill~=2)
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.kill:  ' num2str(VF.kill)]),return
end
if (VF.skip_pole~=0) && (VF.skip_pole)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.skip_pole:  ' num2str(VF.skip_pole)]),return
end
if (VF.skip_res~=0) && (VF.skip_res)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.skip_res:  ' num2str(VF.skip_res)]),return
end
if (VF.use_normal~=0) && (VF.use_normal)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.use_normal:  ' num2str(VF.use_normal)]),return
end
if (VF.use_sparse~=0) && (VF.use_sparse)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.use_sparse:  ' num2str(VF.use_sparse)]),return
end
if (VF.cmplx_ss~=0) && (VF.cmplx_ss)~=1
  disp(['    ERROR in vectfit2.m: ==> Illegal value for VF.cmplx_ss:  ' num2str(VF.cmplx_ss)]),return
end

rmserr=[];%SERC=[];
[a,b]=size(s);
if a<b, s=s.'; end 

% Some sanity checks on dimension of input arrays:
if length(s)~=length(f(1,:))
  disp('Error in vectfit2.m!!! ==> Second dimension of f does not match length of s.'); 
  return;
end  
if length(s)~=length(weight(1,:))
  disp('Error in vectfit2.m!!! ==> Second dimension of weight does not match length of s.');
  return;
end
if length(weight(:,1))~=1
  if length(weight(:,1))~=length(f(:,1))
    disp('Error in vectfit2.m!!! ==> First dimension of weight is neither 1 nor matches first dimension of f.');  
    return;
  end
end  
    


  %set(0,'DefaultLineLineWidth',0.5) ; set(0,'DefaultLineMarkerSize',4) ;
  %clear b; clear C; 
  LAMBD=diag(poles);
  Ns=length(s); N=length(LAMBD); Nc=length(f(:,1));
  B=ones(N,1); %I=diag(ones(1,N));                
  SERA=poles;SERC=zeros(Nc,N);SERD=zeros(Nc,1);SERE=zeros(Nc,1);
  roetter=poles;
  fit=zeros(Nc,Ns);

  weight=weight.';
  if length(weight(1,:))==1
    common_weight=1;
  elseif length(weight(1,:))==Nc
    common_weight=0;
  else
    disp('ERROR in vectfit_new_sigma7.m: Invalid size of array weight')
    return
  end

  if VF.asymp==1
    offs=0; 
  elseif VF.asymp==2
    offs=1;  
  else
    offs=2; 
  end
  
  
  
%=========================================================================
%=========================================================================
%  POLE IDENTIFICATION:
%=========================================================================
%=========================================================================
 
  
if VF.skip_pole~=1  
  

if VF.use_sparse==1
  SS=zeros(Nc*Ns*(N+offs) + Nc*Ns*(N+1) +(N+1),1);
  II=SS;
  JJ=SS;
  Escale=sparse(zeros(1,Nc*(N+offs)+N+1));
else
  A=zeros(Ns*Nc+1,(N+offs)*Nc+N+1); %b=zeros(Ns*Nc+1,1);
  Escale=zeros(1,Nc*(N+offs)+N+1);
end %if use_sparse==1

%=======================================================
% Finding out which starting poles are complex :
%=======================================================
cindex=zeros(1,N);
for m=1:N 
  if imag(LAMBD(m,m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end

%=======================================================
% Building system - matrix :
%=======================================================
 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N);
  for m=1:N
    if cindex(m)==0      %real pole
      Dk(:,m)=1./(s-LAMBD(m,m));
    elseif cindex(m)==1  %complex pole, 1st part
      Dk(:,m)  =1./(s-LAMBD(m,m)) + 1./(s-LAMBD(m,m)');
      Dk(:,m+1)=i./(s-LAMBD(m,m)) - i./(s-LAMBD(m,m)');
    end
  end      
 if VF.asymp==1 || VF.asymp==2    
   Dk(:,N+1)=1; 
 elseif VF.asymp==3   
   Dk(:,N+1)=1; 
   Dk(:,N+2)=s;
 end
    
%Scaling for last row of LS-problem (pole identification)
scale=0;
for m=1:Nc
  if length(weight(1,:))==1
    scale=scale+(norm(weight(:,1).*f(m,:).'))^2;
  else
    scale=scale+(norm(weight(:,m).*f(m,:).'))^2;
  end  
end
scale=sqrt(scale)/Ns; 

 
if VF.use_sparse==1  
    
  %Filling in the block diagonal part of A : 
  koko1=linspace(1,Ns,Ns).';
  for n=1:Nc
      
    if common_weight==1
      weig=weight;
    else
      weig=weight(:,n);
    end  
      
    indblokki=(n-1)*Ns*(N+offs); %offset for hver "blokk" (2Ns*(N+offs))
    indblokkj=(n-1)*(N+offs);  
    dum1=Ns*(n-1);
    for m=1:N+offs
      Dkm=weig.*Dk(:,m); 
      ind1=indblokki+(m-1)*Ns;
      II(ind1+1:ind1+Ns)=dum1+koko1;
      JJ(ind1+1:ind1+Ns)=indblokkj+m;
      SS(ind1+1:ind1+Ns)=Dkm;
    end
  %end  

 %Filling in the right blocks (N+1 columns) :
  indblokkleft=Ns*Nc*(N+offs);
  %for n=1:Nc
      
    indblokki=(n-1)*Ns*(N+1);
    indblokkj=Nc*(N+offs);
    dum1=Ns*(n-1);
    for m=1:N+1
      ind1=indblokki+(m-1)*Ns;
      ind11=ind1+indblokkleft;
      II(ind11+1:ind11+Ns)=dum1+koko1;
      JJ(ind11+1:ind11+Ns)=indblokkj+m;
      SS(ind11+1:ind11+Ns)=-weig.*Dk(:,m).*f(n,:).';
    end  
  end       
  
  %Adding the last row (integral criterion for sigma):
  %!!scale=norm(f) %/Ns;

  ind = Ns*Nc*(N+offs) +Ns*Nc*(N+1);
  for mm=1:N+1  
    II(ind+mm)=Nc*Ns+1;
    JJ(ind+mm)=Nc*(N+offs)+mm;      
    SS(ind+mm)=real(scale*sum(Dk(:,mm)));
  end 

  %Right side:
  b=sparse(2*Nc*Ns+2,1);
  b(Nc*Ns+1)=Ns*scale;
    
  SSSS=sparse(II,JJ,SS,Ns*Nc+1,Nc*(N+offs)+N+1); 
  SSSS=[real(SSSS); imag(SSSS)];
  if VF.relax==1
    for m=1:Nc*(N+offs)+N+1
      Escale(m)=norm(SSSS(:,m),2);
      SSSS(:,m)=SSSS(:,m)./Escale(m);  
    end 
    if VF.use_normal==1
      x=(SSSS.'*SSSS)\(SSSS.'*b);    
    else    
      x=SSSS\b;   
    end
    if VF.relax==1
      x=x./Escale.';
    end  
  else
    x=0;
  end
    

  %*****************October 20, 2005*******************
  if VF.relax==0 | abs(x(end))<TOLlow | abs(x(end))>TOLhigh
    if VF.relax==0
      Dnew=1; 
    else
      if x(end)==0
        Dnew=1;
      elseif abs(x(end))<TOLlow
        Dnew=sign(x(end))*TOLlow
      elseif abs(x(end))>TOLhigh 
        Dnew=sign(x(end))*TOLhigh
      end
      for col=1:length(SSSS(1,:));
        SSSS(:,col)=SSSS(:,col).*Escale(col); %removing previous scaling
      end      
    end  

    ind=length(SSSS(:,1))/2; %index to additional row related to relaxation
    SSSS(ind,:)=[];  
    b=-Dnew*SSSS(:,end);  %new right side 
    SSSS(:,end)=[]; 
    Escale(end)=[];
    for col=1:length(SSSS(1,:));
      Escale(col)=norm(SSSS(:,col),2); %Euclidian norm
      SSSS(:,col)=SSSS(:,col)./Escale(col);
    end
    if VF.use_normal==1      
      x=(SSSS.'*SSSS)\(SSSS.'*b);    
    else
      x=SSSS\b;   
    end
    x=x./Escale.';  
    x=[x;Dnew];
  end %if abs(x(end))<TOL
  %************************************ 

else %if use_sparse==1 

  %scale=norm(f);%/Ns; %Scaling for sigma in LS problem 
  for n=1:Nc
    
    if common_weight==1
      weig=weight;
    else
      weig=weight(:,n);
    end 

    inda=(n-1)*Ns+1; indb=n*Ns;
    ind1=(n-1)*(N+offs)+1; 
    for m=1:N+offs %Diagonal blocks
      A(inda:indb,ind1-1+m)=weig.*Dk(1:Ns,m); 
    end  
    for m=1:N+1 %Rightmost blocks
      A(inda:indb,Nc*(N+offs)+m)=-weig.*(Dk(1:Ns,m).*f(n,1:Ns).');
    end          
    
    
  end %for n=1:Nc 
  
  %Integral criterion for sigma:
  offset=Nc*(N+offs);
  for mm=1:N+1  
    A(Nc*Ns+1,offset+mm)=real(scale*sum(Dk(:,mm)));
  end   
  b=zeros(Nc*Ns+1,1);
  b(Nc*Ns+1)=Ns*scale; 

   xb =real(b); xxb =imag(b); 

  Nrow=Ns*Nc+1; 
  %Ncol=(N+offs)*Nc+N+1;

  b=zeros(2*Nrow,1);
  %A(Nrow+1:2*Nrow,1:Ncol) = imag(A(1:Nrow,1:Ncol));     
  %A(1:Nrow,1:Ncol) = real(A(1:Nrow,1:Ncol));
  A=[real(A);imag(A)];

  b(1:Nrow,1)           = xb;            
  b(Nrow+1:2*Nrow,1)    = xxb;  
  %clear xb; clear xxb;

  if VF.relax==1
    for col=1:length(A(1,:));
      Escale(col)=norm(A(:,col),2); %Euclidian norm
      A(:,col)=A(:,col)./Escale(col);
    end
    if VF.use_normal==1      
      x=(A.'*A)\(A.'*b); 
    else
      x=A\b;       
    end
    x=x./Escale.';
  else
    x=0;  
  end

  %*****************October 20, 2005*******************
  if VF.relax==0 | abs(x(end))<TOLlow | abs(x(end))>TOLhigh
    if VF.relax==0
      Dnew=1; 
    else
      if x(end)==0
        Dnew=1;      
      elseif abs(x(end))<TOLlow
        Dnew=sign(x(end))*TOLlow
      elseif abs(x(end))>TOLhigh 
        Dnew=sign(x(end))*TOLhigh   
      end
      for col=1:length(A(1,:));
        A(:,col)=A(:,col).*Escale(col); %removing previous scaling
      end        
    end
  
    ind=length(A(:,1))/2; %index to additional row related to relaxation
    A(ind,:)=[];  
    b=-Dnew*A(:,end);  %new right side 
    A(:,end)=[]; 
    Escale(end)=[];
    for col=1:length(A(1,:));
      Escale(col)=norm(A(:,col),2); %Euclidian norm
      A(:,col)=A(:,col)./Escale(col);
    end
    if VF.use_normal==1      
      x=(A.'*A)\(A.'*b); 
    else
      x=A\b;       
    end
    x=x./Escale.';
    x=[x;Dnew];
  end %if abs(x(end))<TOL
  
 
  %************************************  
  
  
  
end %if use_sparse==1 






%clear A;
if VF.asymp==1
  for j=1:Nc
    C(:,j)=x((j-1)*(N)+1:j*(N));
  end
  C(:,Nc+1)=x(Nc*(N)+1:Nc*(N)+N);
  C=C.';
elseif VF.asymp==2
  for j=1:Nc
    C(:,j)=x((j-1)*(N+1)+1:j*(N+1)-1);
    SERD(j,1)=x(j*(N+1));
  end
  C(:,Nc+1)=x(Nc*(N+1)+1:Nc*(N+1)+N);
  C=C.';
else
  for j=1:Nc
    C(:,j)=x((j-1)*(N+2)+1:j*(N+2)-2);
    SERD(j,1)=x(j*(N+2)-2+1);
    SERE(j,1) =x(j*(N+2)-2+2);
  end
  C(:,Nc+1)=x(Nc*(N+2)+1:Nc*(N+2)+N);
  C=C.';
end
D=x(end); %NEW!!

%We now change back to make C complex : 
% **************
for m=1:N
  if cindex(m)==1
    for n=1:Nc+1
      r1=C(n,m); r2=C(n,m+1);
      C(n,m)=r1+i*r2; C(n,m+1)=r1-i*r2;
    end
  end
end
% **************  

if VF.spy1==1
  Dk=zeros(Ns,N); 
  for m=1:N
    Dk(:,m)=1./(s-LAMBD(m,m));
  end 
  for n=1:Nc
    RES1(:,n)=Dk*C(n,:).';
    if VF.asymp==2
      RES1(:,n)=RES1(:,n)+SERD(n); 
    elseif VF.asymp==3
      RES1(:,n)=RES1(:,n)+SERD(n)+s.*SERE(n); % (sigma*f)rat 
    end
    RES2(:,n)=(D+Dk*C(Nc+1,:).').*f(n,:).'; % fnum*(sigma)rat
  end 
  RES3(:,1)=D+Dk*C(Nc+1,:).'; % (sigma)rat
  freq=s./(2*pi*i); 
  if VF.logx==1  
    if VF.logy==1
      figure(3),
      h1=loglog(freq,abs(RES1'),'c'); xlim([freq(1) freq(Ns)]);hold on  %sigma*f
      h2=loglog(freq,abs(RES3.'),'r');        %sigma
      h3=loglog(freq,abs(RES2.'-RES1.'),'g');hold off; 
    else %logy=0
      figure(3),
      h1=semilogx(freq,abs(RES1'),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogx(freq,abs(RES3.'),'r');
      h3=semilogx(freq,abs(RES2.'-RES1.'),'g');hold off; pause(0.1);
    end
  else %logx=0
    if VF.logy==1
      figure(3),
      h1=semilogy(freq,abs(RES1'),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogy(freq,abs(RES3.'),'r');
      h3=semilogy(freq,abs(RES2.'-RES1.'),'g');hold off; 
    else %logy=0
      figure(3),
      h1=plot(s./(2*pi*i),abs(RES1'),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=plot(s./(2*pi*i),abs(RES3.'),'r');
      h3=plot(s./(2*pi*i),abs(RES2.'-RES1.'),'g');hold off; 
    end
  end
  figure(3),xlabel('Frequency [Hz]'); ylabel('Magnitude'); 
  figure(3), title('Sigma and sigma*f')
  if VF.legend==1
    legend([h1(1) h2(1) h3(1)],'(sigma*f)_{rat}','sigma','Deviation: (sigma*f)_{rat} - sigma*f');
  end  
  drawnow;
end

%=============================================================================
% We now calculate the zeros for sigma :
%=============================================================================
%oldLAMBD=LAMBD;oldB=B;oldC=C;
m=0;
for n=1:N
  m=m+1;
  if m<N  
    if( abs(LAMBD(m,m))>abs(real(LAMBD(m,m))) ) %complex number?
      LAMBD(m+1,m)=-imag(LAMBD(m,m)); LAMBD(m,m+1)=imag(LAMBD(m,m));
      LAMBD(m,m)=real(LAMBD(m,m));LAMBD(m+1,m+1)=LAMBD(m,m);
      B(m,1)=2; B(m+1,1)=0; 
      koko=C(Nc+1,m); C(Nc+1,m)=real(koko); C(Nc+1,m+1)=imag(koko);
      m=m+1;
    end
  end
end
ZER=LAMBD-B*C(Nc+1,:)/D;
roetter=eig(ZER).'; 
unstables=real(roetter)>0;  
if VF.kill==1
  roetter(unstables)=[]; %Deleting unstable zeros (they become poles later...)
elseif VF.kill==2
  roetter(unstables)=roetter(unstables)-2*real(roetter(unstables)); %Forcing unstable poles to be stable...
end
roetter=sort(roetter); N=length(roetter);


%=============================================
%Sorterer polene s.a. de reelle kommer først:
for n=1:N
  for m=n+1:N
    if imag(roetter(m))==0 && imag(roetter(n))~=0
      trans=roetter(n); roetter(n)=roetter(m); roetter(m)=trans;
    end
  end
end
N1=0;
for m=1:N
  if imag(roetter(m))==0, N1=m; end
end
if N1<N, roetter(N1+1:N)=sort(roetter(N1+1:N)); end   % N1: n.o. real poles
%N2=N-N1;                                             % N2: n.o. imag.poles

roetter=roetter-2*i*imag(roetter); %10.11.97 !!!
SERA=roetter.';


end %if skip_pole~=1  

%=========================================================================
%=========================================================================
%  RESIDUE IDENTIFICATION:
%=========================================================================
%=========================================================================

if VF.skip_res~=1

%=============================================================================
% We now calculate SER for f, using the modified zeros of sigma as new poles :
%========================================================================================

%clear LAMBD A A1 xA1 xxA1 A2 xA2 xxA2 b xb xxb C RES1 RES2;

LAMBD=roetter; 

%B=ones(N,1); 

% Finding out which poles are complex :
cindex=zeros(1,N);
for m=1:N 
  if imag(LAMBD(m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end


%===============================================================================
% We now calculate the SER for f (new fitting), using the above calculated 
% zeros as known poles :
%=============================================================================== 
if VF.asymp==1
  A=zeros(2*Ns,N); BB=zeros(2*Ns,Nc);
elseif VF.asymp==2
  A=zeros(2*Ns,N+1); BB=zeros(2*Ns,Nc);
else
  A=zeros(2*Ns,N+2); BB=zeros(2*Ns,Nc);
end

 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N); 
    for m=1:N
      if cindex(m)==0      %real pole
        Dk(:,m)=1./(s-LAMBD(m));
      elseif cindex(m)==1  %complex pole, 1st part
        Dk(:,m)  =1./(s-LAMBD(m)) + 1./(s-LAMBD(m)');
        Dk(:,m+1)=i./(s-LAMBD(m)) - i./(s-LAMBD(m)');
      end
    end 


if common_weight==1

 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N); 
    for m=1:N
      if cindex(m)==0      %real pole
        Dk(:,m)=weight./(s-LAMBD(m));
      elseif cindex(m)==1  %complex pole, 1st part
        Dk(:,m)  =weight./(s-LAMBD(m)) + weight./(s-LAMBD(m)');
        Dk(:,m+1)=i.*weight./(s-LAMBD(m)) - i.*weight./(s-LAMBD(m)');
      end
    end 

  if VF.asymp==1
    A(1:Ns,1:N)=Dk; 
  elseif VF.asymp==2
    A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=weight;
  else
    A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=weight; A(1:Ns,N+2)=weight.*s;
  end
  for m=1:Nc
    BB(1:Ns,m)=weight.*f(m,:).';
  end
  A(Ns+1:2*Ns,:)=imag(A(1:Ns,:));
  A(1:Ns,:)=real(A(1:Ns,:));
  BB(Ns+1:2*Ns,:)=imag(BB(1:Ns,:));
  BB(1:Ns,:)=real(BB(1:Ns,:));

  if VF.asymp==2
    A(1:Ns,N+1)=A(1:Ns,N+1);             
  elseif VF.asymp==3
    A(1:Ns,N+1)=A(1:Ns,N+1);             
    A(Ns+1:2*Ns,N+2)=A(Ns+1:2*Ns,N+2);  
  end

  %clear Escale;
  Escale=zeros(1,length(A(1,:)));
  for col=1:length(A(1,:));
    Escale(col)=norm(A(:,col),2);
    A(:,col)=A(:,col)./Escale(col);
  end
  X=A\BB;
  for n=1:Nc
    X(:,n)=X(:,n)./Escale.';
  end

  %clear A;
  X=X.';
  C=X(:,1:N); 
  if VF.asymp==2
    SERD=X(:,N+1);
  elseif VF.asymp==3
    SERE=X(:,N+2); 
    SERD=X(:,N+1);
  end

else %if common_weight==1
    
  SERD=zeros(Nc,1); 
  SERE=zeros(Nc,1); 
  C=zeros(Nc,N);
  for n=1:Nc

    if VF.asymp==1
      A(1:Ns,1:N)=Dk; 
    elseif VF.asymp==2
      A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=1;
    else
      A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=1; A(1:Ns,N+2)=s;
    end

    for m=1:length(A(1,:))
      A(1:Ns,m)=weight(:,n).*A(1:Ns,m);
    end
 
    BB=weight(:,n).*f(n,:).';
    A(Ns+1:2*Ns,:)=imag(A(1:Ns,:));
    A(1:Ns,:)=real(A(1:Ns,:));
    BB(Ns+1:2*Ns)=imag(BB(1:Ns));
    BB(1:Ns)=real(BB(1:Ns));

    if VF.asymp==2
      A(1:Ns,N+1)=A(1:Ns,N+1);             
    elseif VF.asymp==3
      A(1:Ns,N+1)=A(1:Ns,N+1);             
      A(Ns+1:2*Ns,N+2)=A(Ns+1:2*Ns,N+2);  
    end

    %clear Escale;
    Escale=zeros(1,length(A(1,:)));
    for col=1:length(A(1,:));
      Escale(col)=norm(A(:,col),2);
      A(:,col)=A(:,col)./Escale(col);
    end
    x=A\BB;
    x=x./Escale.';

    %clear A;
    C(n,1:N)=x(1:N).'; 

    if VF.asymp==2
      SERD(n)=x(N+1);
    elseif VF.asymp==3
      SERE(n)=x(N+2); 
      SERD(n)=x(N+1);
    end


  end %for n=1:Nc

    
end %if common_weight==1    
  

%=========================================================================

%We now change back to make C complex. 
for m=1:N
  if cindex(m)==1
    for n=1:Nc
      r1=C(n,m); r2=C(n,m+1);
      C(n,m)=r1+i*r2; C(n,m+1)=r1-i*r2;
   end
  end
end
% **************  


B=ones(N,1);

%====================================================

SERA  = LAMBD;
SERB  = B;
SERC  = C;



  Dk=zeros(Ns,N); 
  for m=1:N
    Dk(:,m)=1./(s-SERA(m));
  end 
  for n=1:Nc
    fit(n,:)=(Dk*SERC(n,:).').';
    if VF.asymp==2
      fit(n,:)=fit(n,:)+SERD(n);
    elseif VF.asymp==3
      fit(n,:)=fit(n,:)+SERD(n)+s.'.*SERE(n);
    end
  end 
  
  fit=fit.'; f=f.';
  diff=fit-f; rmserr=sqrt(sum(sum(abs(diff.^2))))/sqrt(Nc*Ns);

if VF.spy2==1
  freq=s./(2*pi*i);
  if VF.logx==1     
    if VF.logy==1
      figure(1),
      h1=loglog(freq,abs(f),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=loglog(freq,abs(fit),'m--'); hold off
      if VF.errplot== 1
        hold on,h3=loglog(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1),
      h1=semilogx(freq,abs(f),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogx(freq,abs(fit),'m--'); hold off
      if errplot== 1
        hold on,h3=semilogx(freq,abs(f-fit),'g');hold off; 
      end
    end
    if VF.phaseplot==1
      figure(2),
      h4=semilogx(freq,180*unwrap(angle(f))/pi,'c'); xlim([freq(1) freq(Ns)]);hold on
      h5=semilogx(freq,180*unwrap(angle(fit))/pi,'m--');hold off
    end  
  else %logx=0
    if VF.logy==1
      figure(1),
      h1=semilogy(freq,abs(f),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogy(freq,abs(fit),'m--'); hold off
      if VF.errplot== 1
        hold on,h3=semilogy(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1), 
      h1=plot(freq,abs(f),'c'); xlim([freq(1) freq(Ns)]);hold on
      h2=plot(freq,abs(fit),'m--'); hold off
      if VF.errplot== 1
        hold on,h3=plot(freq,abs(f-fit),'g');hold off; 
      end
    end
    if VF.phaseplot==1    
      figure(2),
      h4=plot(freq,180*unwrap(angle(f))/pi,'c'); xlim([freq(1) freq(Ns)]);hold on
      h5=plot(freq,180*unwrap(angle(fit))/pi,'m--');hold off
    end  
  end %logy=0
  figure(1),
  xlabel('Frequency [Hz]'); ylabel('Magnitude [p.u.]');
  title('Approximation of f');
  if VF.legend==1
    if VF.errplot==1
      legend([h1(1) h2(1) h3(1)],'Original','VF','Deviation');
    else
      legend([h1(1) h2(1)],'Original','VF');       
    end  
  end
  if VF.phaseplot==1
    figure(2),
    xlabel('Frequency [Hz]'); ylabel('Phase angle [deg]');
    title('Approximation of f');
    if VF.legend==1
      legend([h4(1) h5(1)],'Original','VF'); 
    end
  end
  drawnow;
end
fit=fit.';


end %if skip_res~=1

A=SERA;
poles=A;
if VF.skip_res~=1
  B=SERB; C=SERC; D=SERD; E=SERE;
else
  B=ones(N,1); C=zeros(Nc,N); D=zeros(Nc,Nc); E=zeros(Nc,Nc); rmserr=0;  
end;    


%========================================
% Convert into real state-space model
%========================================
if VF.cmplx_ss~=1 

  A=diag(sparse(A));
  
  cindex=zeros(1,N);
  for m=1:N 
    if imag(A(m,m))~=0  
      if m==1 
        cindex(m)=1;
      else
        if cindex(m-1)==0 || cindex(m-1)==2
          cindex(m)=1; cindex(m+1)=2; 
        else
          cindex(m)=2;
        end
      end 
    end
  end
  
  n=0;
  for m=1:N
    n=n+1;
    if cindex(m)==1
      a=A(n,n); a1=real(a); a2=imag(a);
      c=C(:,n); c1=real(c); c2=imag(c);
      b=B(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
   
      A(n:n+1,n:n+1)=Ablock;
      C(:,n)=c1;
      C(:,n+1)=c2;
      B(n,:)=b1;
      B(n+1,:)=b2;
    end
  end

else
  A=sparse(diag(A));    % A is complex, make it diagonal
end %if cmplx_ss~=1    
    
SER.A=A; SER.B=B; SER.C=C; SER.D=D; SER.E=E;


