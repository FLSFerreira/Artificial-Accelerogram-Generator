% ARTIFICIAL ACCELEROGRAM GENERATOR CODE (AAG) WRITTEN BY FERREIRA et al. (2020)
function [art_sign,art_spec,Sa]=AAG(soil,type,importance,agr)
% INPUT (Eurocode 8 section 3.2.2.2)
% soil=1 to 5 equivalent to soil A to E
% Earthquake type=1 or type=2
% importance = 0.8, 1, 1.2, 1.4 respectively for importance class I to IV
% agr is the base peak acceleration reference
% OUTPUT
% art_sign is the artificial sign and art_spec is the artificial spectra
load Real_Seismic_Events.mat Real_Seismic_Events
%Real_Seismic_Events 1st line is the time vector,
%2nd line the ElCentro Earthquake (1940), 3th line the Gebze Earthquake (1999)
%4th line the Mexico City Earthquake (1985)
%USER DEFINED PARAMETERS
T=0.01:0.02:4; %Period range to define spectra
z=0.05; %Intended spectra damping value z=5% for EC8
N_it=8; %number of iterations
calculate_events=1:3; % if the user wants to use only some events
% all earthquakes are trimmed at 40.96s with a dt=0.005s (8192 points=2^13)
number_points=length(Real_Seismic_Events(1,:)); %MUST BE a power of 2
% Elastic spectra definition
[Sa]=EC8elasticspectra(soil,type,importance,z,agr);
[Sa_T]=[T; interp1(Sa(1,:),Sa(2,:),T)];
% Cicle to generate 3 artificial acelerograms
art_sign=zeros(3,number_points);
art_spec=zeros(3,length(T));
aux_vE=zeros(4,number_points);
time=Real_Seismic_Events(1,:);
aux_vE(1,:)=time;
for jj=calculate_events
aux_vE(2,:)=Real_Seismic_Events(jj+1,:);
[Final_accelerogram,S]=F_AAG(aux_vE,N_it,Sa_T,T,z);
art_sign(jj,:)=Final_accelerogram(2,:); %The final acelerogram is saved
art_spec(jj,:)=S(end,:);
end
% Artificial signal and spectra plot
for j=calculate_events %Plot
figure(j)
plot(Sa_T(1,:),Sa_T(2,:),'b')
hold on
plot(T,art_spec(j,:),'r')
figure(100+j)
plot(time,art_sign(j,:))
end
art_sign=[time;art_sign];  
art_spec=[T;art_spec];
end
%----AUXILIARY FUNCTIONS----
% GENERATE EUROCODE 8 HORIZONTAL ELASTIC SPECTRA 
function [S_elastic]=EC8elasticspectra(soil,type,importance,z,agr)
division=100;
Ttrim=10; %period above which the spectra is trimed
if type==1
          %  S   TB   TC  TD %  S    TB   TC  TD %  S    TB   TC  TD
aux_Matrix=[1,   0.15, 0.4, 2; 1.2, 0.15, 0.5, 2; 1.15, 0.2, 0.6, 2;
            1.35, 0.2, 0.8, 2; 1.4, 0.15, 0.5, 2];
elseif type==2
               %S     TB     TC    TD
aux_Matrix=[1, 0.05, 0.25, 1.2; 1.35, 0.05, 0.25, 1.2; 1.5, 0.1, 0.25, 1.2;
            1.8, 0.1, 0.3, 1.2; 1.6,  0.05, 0.25, 1.2];
end
S=aux_Matrix(soil,1); TB=aux_Matrix(soil,2); TC=aux_Matrix(soil,3);
TD=aux_Matrix(soil,4); ag=importance*agr;
damping_modif=max(sqrt(10/(5+100*z)),0.55);
TC_TD=(TC):(TD-TC)/division:TD;
aux_TC_TD=ag*S*damping_modif*2.5*TC./TC_TD;
TD_Ttrim=(TD+(4-TD)/division):(Ttrim-TD)/division:Ttrim;
aux_TD_Ttrim=ag*S*damping_modif*2.5*TC*TD./(TD_Ttrim.^2);
% The final result has single points in the linear region 
S_elastic=[0,    TB,                    TC_TD,     TD_Ttrim;
           S*ag, S*ag*damping_modif*2.5, aux_TC_TD, aux_TD_Ttrim];
end       
% TRANSFORM THE ORIGINAL EVENT TO A RESPONSE SPECTRA MATCHED EVENT    
function [Final_accelerogram,Sa]=F_AAG(original_accel,N_iter,Sa_elastic,T,z)
% INPUT
% original_accel is the original accelerogram to be transformed:
% original_accel = [time; acceleration; velocity; displacement]
% N_iter is the number of iterations of the algorithm 
% Sa_elastic is the elastic spectra to be matched
% T is the period vector
% z is the structural damping
% ag is the objective peak ground acceleration (PGA)
%OUTPUT
% Final_accelerogram is the final artificial accelerogram in the format:
% Final accelerogram = [time; acceleration; velocity; displacement]
% Sa is the acceleration response spectra in all iterations
% S(N_iter,:) is the Final_accelerogram response spectra
ag=Sa_elastic(2,1); %Peak Ground Acceleration
freq_spectra=1./T;
Sa_design_interp=interp1(Sa_elastic(1,:),Sa_elastic(2,:),T);
% Scaling of original accelerogram
sfac=0.85; %Scaling in the first iteration so that the PGA=sfac*ag 
original_accel(2,:)=original_accel(2,:)*ag/max(abs(original_accel(2,:)))*sfac;
time=original_accel(1,:);
y=original_accel(2,:);
dt=original_accel(1,2)-original_accel(1,1);
Sa=zeros(N_iter,length(T));
%Auxiliar variables (!!number_points)
number_points=length(original_accel(1,:)); %MUST BE a power of 2
aux_var_1=[0, 0.05, 0.95,1;
           0,    1,    1,0];
aux_var_2=interp1(aux_var_1(1,:),aux_var_1(2,:),0:1/(number_points-1):1); 
P_MODIF=zeros(N_iter,length(T)); %accelerogram modification factors 
SCALE_factor=zeros(N_iter,1); %accelerogram scale factors 
%original_accel response spectra
S_aux=ARS([time;y],z,T); 
Sa(1,:)=S_aux(4,:);
P_MODIF(1,:)=Sa_design_interp./Sa(1,:);
SCALE_factor(1,:)=sum(P_MODIF(1,:))/length(P_MODIF(1,:));
P_MODIF(1,:)=P_MODIF(1,:)/SCALE_factor(1);
% Accelerogram Fourier transform (FFT)
L=length(y);
Fs=1/dt;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs*linspace(0,1,NFFT);
% Accelerogram FFT modification to match the objective spectra
MODIF_Y=zeros(N_iter-1,NFFT);
Y_modif=zeros(N_iter-1,NFFT);
V_modif=zeros(N_iter-1,NFFT);
D_modif=zeros(N_iter-1,NFFT);
y_aux=zeros(N_iter,L);
v_aux=zeros(N_iter,L);
d_aux=zeros(N_iter,L);
y_aux(1,1:L)=original_accel(2,:);
v_aux(1,1:L)=original_accel(3,:);
d_aux(1,1:L)=original_accel(4,:);
for j=2:N_iter
P_MODIF_aux=P_MODIF(j-1,:).*SCALE_factor(j-1);         
for jj=1:(NFFT/2)
    if f(jj)>=freq_spectra(end) && f(jj)<=freq_spectra(1)
        MODIF_Y(j-1,jj)=interp1(freq_spectra,P_MODIF_aux,f(jj));
    else
        MODIF_Y(j-1,jj)=1;
    end    
end
MODIF_Y(j-1,NFFT:-1:(NFFT/2+1))=MODIF_Y(j-1,1:(NFFT/2)); %Aliasing keeping
if j==2 %1st iteration
Y_modif(j-1,:)=Y.*MODIF_Y(j-1,:);
V_modif(j-1,:)=Y_modif(j-1,:)./(i*2*pi*f);
D_modif(j-1,:)=-Y_modif(j-1,:)./(2*pi*f).^2;
else 
Y_modif(j-1,:)=Y_modif(j-2,:).*MODIF_Y(j-1,:);
V_modif(j-1,:)=Y_modif(j-1,:)./(i*2*pi*f);
D_modif(j-1,:)=-Y_modif(j-1,:)./(2*pi*f).^2;
end
D_modif(j-1,1)=0; % first frequency = 0 is not intended to be considered 
V_modif(j-1,1)=0;
% Modified accelerogram determination
y_recomposto_modif_aux = real(ifft(Y_modif(j-1,:))*L);
y_aux(j,:)=y_recomposto_modif_aux(1:L).*aux_var_2;
velo_recomposto_modif_aux = real(ifft(V_modif(j-1,:))*L);
v_aux(j,:)=velo_recomposto_modif_aux(1:L);
desl_recomposto_modif_aux = real(ifft(D_modif(j-1,:))*L);
d_aux(j,:)=desl_recomposto_modif_aux(1:L)-desl_recomposto_modif_aux(1);
% Modified accelerogram response spectra
S_aux=ARS([time;y_aux(j,1:L)],z,T);
Sa(j,:)=S_aux(4,:);
P_MODIF(j,:)=Sa_design_interp./Sa(j,:);
SCALE_factor(j,:)=sum(P_MODIF(j,:))/length(P_MODIF(j,:)); 
P_MODIF(j,:)=P_MODIF(j,:)/SCALE_factor(j);
end
% Amp=[abs(Y);abs(Y_modif)]; %useful variables
% FFT_iter=[P;Amp];
Final_accelerogram=[time;y_aux(N_iter,1:L);v_aux(N_iter,1:L);d_aux(N_iter,1:L)];
end
% DETERMINE THE ACCELEROGRAM RESPONSE SPECTRA (ARS)
function [S]=ARS(accelerogram,z,T)
% Function to determine the linear elastic spectra of a SDOF oscillator 
% INPUT
% T is the oscillator periods vector (s)
% z is the oscillator damping
% accelerogram is the input excitation: [time; acceleration]
% OUTPUT
% S is the response spectra: in the form: 
% S= [period; (displacement, velocity and acceleration response spectra)]
S=zeros(4,length(T));
S(1,:)=T;
w=1./T*2*pi;
%Total_time=accelerogram(1,end);
dt=accelerogram(1,2)-accelerogram(1,1);
time=accelerogram(1,:);
Forces=interp1(accelerogram(1,:),accelerogram(2,:),time);
dForces_dt=Forces(2:end)-Forces(1:(end-1));
for jj=1:length(w); %parfor using parallel computing toolbox
% State Space Differential Equations Matrixes
A=[ 0         1;
    -w(jj)^2, -2*z*w(jj)];
B=[ 0; 1];
% Analytical Time Step Solver of State Space 
Kb=real(exp(1)^(A*dt));
Kf=A\(Kb-eye(2));
KfB=Kf*B;
KdfB=A\(Kf/dt-eye(2))*B;
% Diffenrential Equation Solver
x=zeros(2,length(time)); %1st line -displacement; 2nd line - velocity
acel=zeros(1,length(time)); % acceleration 
aux_acel=-[w(jj)^2,2*z*w(jj)]; 
for j=1:(length(time)-1)   
   x(:,j+1)=Kb*x(:,j)+KfB*Forces(1,j)+KdfB*dForces_dt(1,j);
   acel(j+1)=aux_acel*x(:,j+1); %acceleration
end
%Displacement, velocity and acceleration response spectra
S(2:4,jj)=[max(abs(x(1,:)));max(abs(x(2,:)));max(abs(acel))]; 
end
end
%Arias Intensity
%IA=pi()/(2*9.8)*(sum(art_sign(2:end,:)'.^2)')*art_sign(1,2); 
% =========================================================================
% === This code was written by Ferreira F, Moutinho C, Cunha A and      ===
% === Caetano E.                                                        ===
% === ----------------------------------------------------------------- ===
% === The author can be contacted at              : fferreira@uc.pt     ===
% === ----------------------------------------------------------------- ===
% === The code is intended for research and academic purposes           ===
% === the details can be found in the paper:                            ===
% === Ferreira F. et al, "An Artificial Accelerogram Generator Code     ===
% === Written in Matlab", Engineering Reports. 2020; 2(3):e12129        ===
% === doi: 10.1002/eng2.12129                                           ===
% === ----------------------------------------------------------------- ===
% === The code and data can be downloaded from :                        ===
% === https://github.com/FLSFerreira/Artificial-Accelerogram-Generator  ===
% === and in Supporting Information at :                                ===
% === https://onlinelibrary.wiley.com/doi/10.1002/eng2.12129            ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed, modified and used for academic purposes=
% === The paper should be properly and appropriately referenced         ===
% === The author do not guarantee that the code is free from errors, and ==
% === does not take any responsability by the use of the code           ===


