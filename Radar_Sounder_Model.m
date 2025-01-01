%% Driver Program to Apply Numerical Method to Radar Sounder Model
%% Kirchhoff Approximation
%% Incident beam is Gaussian beam 
%% Incident pules is Gaussian pulse
%% incident angle = or ~ 0
%% Updates: Is_0 is calculated by E*conj(E). Input paremeters are put into a matrix P.
%%                Long equations are divided into shot equations.
%%%% Inputs
%% freq = frequency
%% tau = delayed time
%% d = depth
%% ep2r = dielectric constant of the 2nd layer
%% ep3r = dielectric constant of the 3rd layer
%% sig1 = Top Layer Roughness RMS height
%% L1 = Top Layer Roughness Correlation Length
%% sig2 = Bottom Layer Roughness RMS height
%% L2 = Bottom Layer Roughness Correlation Length
%% plot top and bottom vs time

clc;clear;
c0=299792458; 
gx=0.5550;
gy=9.97e-5;
freq=930e6; % frequency, L
w0=2*pi*freq;
tp = 1.67e-9; % tp is tau_p in equation (17) of (Picardi, 2004).  pulse length = 4*tp
theta0=0;

sig1= 1 /100; % surface variance 
L1= 30 /100; % correlation length
sig2_cm = 1;
sig2= sig2_cm /100; % surface variance 
L2= 20 /100; % correlation length
h0=250000; % height

real_2=3.5;    % real part of ep2r
ep2r=real_2-real_2*0.005i; % layer 2, dry regolith, L
ep3r=3.1500 - 0.00064i ;  % layer 3, ice, L

P1=zeros(1,10); % array of input parameters upper surface
P1(1)=w0;  P1(2)=tp;  P1(3)=gx;  P1(4)=gy;  P1(5)=h0; 
P1(7)=sig1;  P1(8)=L1; P1(9)=c0;  

c2=c0/real(sqrt(ep2r));
P2=zeros(1,10); % array of input parameters lower surface
P2(1)=w0;  P2(2)=tp;  P2(3)=gx;  P2(4)=gy;  P2(5)=h0; 
P2(7)=sig2; P2(8)=L2;  P2(9)=c2;

tau= (0*tp): ( 0.5* tp) : (20*tp); % time difference
d=2; % layer depth

%% Functions defined in this program.

% [Is0, Is1, Is2, Is3, Is_inco, Is]=NM(tau,P)
% apply numerical method to surface scattering 

% [atv,ath,a,delay,Tv12,Tv21,Th12,Th21,Rv12,Rh12,Rv23,Rh23]=RL(freq,d,theta1,ep2r,ep3r,sig)
% calculate coefficients of Regolith Layer
%%
A=pi/(4*gx*gy);
cs=4*pi*h0^4/A; % cross section
cs_dB=dB(cs);

[~,ath,~,delay,~,~,~,~,~,Rh12,~,Rh23]=RL(freq,d,theta0,ep2r,ep3r,sig1);
P1(10)=Rh12;    P2(10)=Rh23;
attn_dB=dB(ath);

%function[Is0, Is1, Is2, Is3, Is_inco, Is]=NM(tau,P)
for J=1:length(tau)

[top0(J),  top1(J), top2(J), top3(J), top_i(J), top(J) ]=NM(tau(J),  P1 ); 
top_dB(J)=dB(top(J))+cs_dB;
top0_dB(J)=dB(top0(J))+cs_dB;
topi_dB(J)=dB(top_i(J))+cs_dB;

[s2_0(J),  s2_1(J), s2_2(J), s2_3(J), s2_i(J),  s2(J)  ]=NM(tau(J)-delay,  P2 ); 
bot_0_dB(J)=dB(s2_0(J))+attn_dB+cs_dB; 
bot_i_dB(J)=dB(s2_i(J))+attn_dB+cs_dB; 
bot_dB(J)=dB(s2(J))+attn_dB+cs_dB;

index=J
end

for n=1:length(tau);
   if top0_dB(n)<-80
      top0_dB(n)=-inf;
   end
   if bot_0_dB(n)<-80
      bot_0_dB(n)=-inf;
   end
end

figure; plot(tau, top0_dB , tau, topi_dB, tau , bot_0_dB ,tau, bot_i_dB );
grid on;

%% Data in dB
function [DB]=dB(x)
DB=10*log10(x);
end

%% coefficients of Regolithh Layer
function [atv,ath,a,delay,Tv12,Tv21,Th12,Th21,Rv12,Rh12,Rv23,Rh23]...
    =RL(freq,d,theta1,ep2r,ep3r,sig)
% a = attenuation (1-way power loss is double wave loss)
% atv, ath = 2-way attenuation, transmissivity included
% delay = 2-way time delay
% T = transmisivity
% R = reflectivity

% parameters
c0=299792458; %speed of light 

% wavenumbers
k0 = 2.0*pi*freq/c0; % wavenumber in free space
k2=k0*sqrt(ep2r); % wave number in 2nd layer
%k3=k0*sqrt(ep3r); % wave number in 3rd layer
rk2=real(k2); ik2=-imag(k2); % k2=k2'-i*k2"
%rk3=real(k3); ik3=-imag(k3); % k3=k3'-i*k3"
%iep2r=-imag(ep2r); % real part-epsilon-2nd medium-relative
a=2*ik2; % 1-way power loss is double wave loss
ad=a*d;

% angles, transmit coefficient, reflection coefficient
c1 = cosd(theta1);
s1 = sind(theta1);
%cf1 = cosd(ph1);
%sf1 = sind(ph1); 

alf = imag(sqrt(ep2r));
be = real(sqrt(ep2r));
p=2*alf*be;

q = be^2 - alf^2 - s1.^2; 
tan2 = 1.414*s1./sqrt(q + sqrt(p.^2 + q.^2));
theta2 = atan(tan2)*180/pi; 
s2 = sind(theta2);

cs2 = sqrt(1 - s2.^2);
tano = 1.414*s1./sqrt(q + sqrt(p.^2 + q.^2)); 
thetao = atan(tano)*180/pi; 
so = sind(thetao);
cso = sqrt(1 - so.^2);

rt = sqrt(ep2r - s1.^2); 
rvi = (ep2r*c1 - rt)/(ep2r*c1 + rt);
rhi = (c1 - rt)/(c1 + rt); 
Rv12 = abs(rvi).^2;
Rh12 = abs(rhi).^2;
rts = sqrt(ep2r - s1.^2); 
rvs = (ep2r*c1 - rts)/(ep2r*c1 + rts); 
rhs = (c1 - rts)/(c1 + rts); 
Rvs = abs(rvs).^2; 
Rhs = abs(rhs).^2;

rtbt = sqrt(ep3r/ep2r - s2.^2); 
rvt = ((ep3r/ep2r)*cs2 - rtbt)/((ep3r/ep2r)*cs2 + rtbt); 
rht = (cs2 - rtbt)/(cs2 + rtbt); 
Rv23 = abs(rvt).^2; 
Rh23 = abs(rht).^2;
%rtbo = sqrt(ep3r/ep2r - so.^2); 
%rvo = ((ep3r/ep2r)* cso - rtbo)/((ep3r/ep2r)* cso + rtbo); 
%rho = (cso - rtbo)/(cso + rtbo); 
%Rvo = abs(rvo)^2; Rho = abs(rho)^2;
Tv12 = 1 - Rv12; Th12 = 1 - Rh12;
Tv21 = 1 - Rvs; Th21 = 1 - Rhs;

atv = exp( -2*(k0-rk2)^2*sig^2 )*Tv12*Tv21*exp(-tau*(1/cso + 1/cs2)); % total attenuation_v
ath = exp( -2*(k0-rk2)^2*sig^2 )*Th12*Th21*exp(-tau*(1/cso + 1/cs2)); % total attenuation_h
v2=2*pi*freq/rk2;
delay=2*d/v2;
end

%% Kirchhoff approximation with Gaussian beam
% Numerical method
% Is0 = 0th order backscattering intensity 
% Is1 = 1st order backscattering intensity
function[Is0, Is1, Is2, Is3, Is_inco, Is]=NM(tau,P)
%% integral of w1, w2 for 0th order
w0=P(1);   tp=P(2);  h0=P(5); 
theta0=P(6); sig=P(7);   c0=P(9);  R=P(10);

ct=cosd(theta0);
R0=h0/ct;
B=R*ct.^2/(16*pi^4*c0^2*R0^4); %  coefficients

wc=w0/(1+2*ct^2*sig^2/(c0^2*tp^2)); % predicted peak of the intergrand
min1=0.3*wc;max1=1.7*wc; % range of integral

sn0=500; % sn number of steps
pieceInt=zeros(1,sn0); 
step1=(max1-min1)/sn0;

 for a=1:sn0
   part1= integrand0( min1+(a-1)*step1 ,tau, P );
   part2= 4* integrand0( min1 + (a-1)*step1+step1/2  ,tau, P );
   part3= integrand0( min1+a*step1, tau, P );
   pieceInt(a) =step1*  (part1+part2+part3) /6 ;  % Simpson integral has 3 parts at each point
 end  
 I0=sum(pieceInt,"all");

 %% integral of w1, w2 for 1st 2nd 3rd order
min1=0.8*wc; max1=1.2*wc; % range of integral
min2=min1;      max2=max1; % range of integral

sn1=500; % sn number of steps
pieceInt1=zeros(sn1,sn1); 
pieceInt2=zeros(sn1,sn1);
pieceInt3=zeros(sn1,sn1);
w1=zeros(1,sn1);
w2=zeros(1,sn1);
step1=(max1-min1)/sn1;
step2=(max2-min2)/sn1;
 for a=1:sn1
  for b=1:sn1
      w1(a)=min1 + (a-1)*step1+step1/2;
      w2(b)=min2 + (b-1)*step2+step2/2 ;
   pieceInt1(a,b) = step1*step2*integrand1( w1(a), w2(b), tau, P)  ;   
   pieceInt2(a,b) = step1*step2*integrand2( w1(a), w2(b), tau, P)  ;
   pieceInt3(a,b) = step1*step2*integrand3( w1(a), w2(b), tau, P)  ;
  end
 end  
%mesh(w1,w2,real(pieceInt1));
I1=sum(pieceInt1,"all");
I2=sum(pieceInt2,"all");
I3=sum(pieceInt3,"all");
I_inco=I1+I2+I3;

 %% results
Is0=B*I0*conj(I0); % 0th order backscattering intensity 

Is1=real(B*I1);

Is2=real(B*I2);

Is3=real(B*I3);

Is_inco=real(B*I_inco);

Is=Is0+Is_inco;

end

%% 0th order integrands of inverse F-transform
% fun0 = 0th order integrand
function[fun0]=integrand0(w1,tau,P)
w0=P(1);   tp=P(2);  gx=P(3);  gy=P(4);  h0=P(5); 
theta0=P(6); sig=P(7);   c0=P(9);  
   
ct=cosd(theta0); st=sind(theta0);

R0=h0/ct;
 
pulse=2*sqrt(pi)*tp*exp(  -tp.^2.*(w1-w0).^2  + 1i.*(w1).*tau );
   
perturb=exp(-2.*ct^2.*sig.^2.*(w1.^2)./c0^2 );
   
ax1=gx^2+1i.*w1/(R0*c0)  + sig^2*sind(2*theta0)^2*w1.^2/(2*R0^2*c0^2); % denominator, x, plus, 0th order
ay1=gy^2+1i.*w1/(R0*c0); % denominator, y, minus, 0th order  
D0=( ax1 .* ay1 ).^0.5; % denominator
   
beta1=-1i*2*st.*w1/c0  + (2*sig^2.*w1.^2*ct*sind(2*theta0))/(c0^2*R0) ;
oblique=exp(  beta1^2/(4*ax1)  );
   
fun0= w1.* pulse * perturb * oblique * pi ./ D0 ; 

end
 
%% 1st order integrands of inverse F-transform
% fun1 = 1st order integrand
function[fun1]=integrand1(w1,w2,tau,P)
 w0=P(1);   tp=P(2);  gx=P(3);  gy=P(4);  h0=P(5); 
 theta0=P(6); sig=P(7);  L=P(8); c0=P(9);  
   ct=cosd(theta0) ;
   R0=h0/ct;
   pulse1=2*sqrt(pi)*tp*exp(  -tp.^2.*(w1-w0).^2+1i.*(w1).*tau );
   pulse2=2*sqrt(pi)*tp*exp( - tp.^2.*(w2-w0).^2 +1i.*(-w2).*tau) ;
   perturb=4*sig^2*(w1*w2/c0^2)*ct^2*exp(    -2.*sig.^2.*ct^2.*(w1.^2+w2.^2)./c0^2     );

   alphax1=gx^2+1i.*w1/(R0*c0)+1/(L^2); 
   alphax2=gx^2-1i.*w2/(R0*c0)+1/(L^2); 
   betax1=-1i*2*sind(theta0).*w1/c0;
   betax2=1i*2*sind(theta0).*w2/c0;
   alphay1=gy^2+1i.*w1/(R0*c0)+1/L^2;
   alphay2=gy^2-1i.*w2/(R0*c0)+1/L^2;
   
   D1=( alphax1*alphax2-1/(L^4) ).^0.5   .*    (alphay1*alphay2-1/(L^4) )^0.5 ;
   oblique=exp( (betax1^2*alphax2+2*betax1*betax2/L^2+betax2^2*alphax1)/( 4*alphax1*alphax2-4/(L^4)) );
 
   fun1=w1*w2*pulse1*pulse2*perturb*oblique*pi^2/D1;
   
 end

%% 2nd order integrands of inverse F-transform
% fun2 = 2nd order integrand
function[fun2]=integrand2(w1,w2,tau,P)
 w0=P(1);   tp=P(2);  gx=P(3);  gy=P(4);  h0=P(5); 
 theta0=P(6); sig=P(7);  L=P(8); c0=P(9);  
   ct=cosd(theta0) ;
   R0=h0/ct;
   pulse1=2*sqrt(pi)*tp*exp(  -tp.^2.*(w1-w0).^2+1i.*(w1).*tau );
   pulse2=2*sqrt(pi)*tp*exp( - tp.^2.*(w2-w0).^2 +1i.*(-w2).*tau) ;
   perturb=8*sig^4*(w1^2*w2^2/c0^4)*ct^4 * exp(    -2.*sig.^2.*ct^2.*(w1.^2+w2.^2)./c0^2     );

   alphax1=gx^2+1i.*w1/(R0*c0)+2/(L^2); 
   alphax2=gx^2-1i.*w2/(R0*c0)+2/(L^2); 
   betax1=-1i*2*sind(theta0).*w1/c0;
   betax2=1i*2*sind(theta0).*w2/c0;
   alphay1=gy^2+1i.*w1/(R0*c0)+2/L^2;
   alphay2=gy^2-1i.*w2/(R0*c0)+2/L^2;
   
   D2=( alphax1*alphax2-4/(L^4) ).^0.5 .* (alphay1*alphay2-4/(L^4) )^0.5 ;
   oblique=exp( (betax1^2*alphax2+4*betax1*betax2/L^2+betax2^2*alphax1)/( 4*alphax1*alphax2-16/(L^4)) );
 
   fun2=w1*w2*pulse1*pulse2*perturb*oblique*pi^2/D2;
end

%% 3rd order integrands of inverse F-transform
% fun3 = 3rd order integrand
function[fun3]=integrand3(w1,w2,tau,P)
 w0=P(1);   tp=P(2);  gx=P(3);  gy=P(4);  h0=P(5); 
 theta0=P(6); sig=P(7);  L=P(8); c0=P(9);  
   ct=cosd(theta0) ;
   R0=h0/ct;
   pulse1=2*sqrt(pi)*tp*exp(  -tp.^2.*(w1-w0).^2+1i.*(w1).*tau );
   pulse2=2*sqrt(pi)*tp*exp( - tp.^2.*(w2-w0).^2 +1i.*(-w2).*tau) ;
   perturb=(32/3)*sig^6*(w1^3*w2^3/c0^6)*ct^6 * exp( -2.*sig.^2.*ct^2.*(w1.^2+w2.^2)./c0^2 );

   alphax1=gx^2+1i.*w1/(R0*c0)+3/(L^2); 
   alphax2=gx^2-1i.*w2/(R0*c0)+3/(L^2); 
   betax1=-1i*2*sind(theta0).*w1/c0;
   betax2=1i*2*sind(theta0).*w2/c0;
   alphay1=gy^2+1i.*w1/(R0*c0)+3/L^2;
   alphay2=gy^2-1i.*w2/(R0*c0)+3/L^2;
   
   D3=( alphax1*alphax2-9/(L^4) ).^0.5 .* (alphay1*alphay2-9/(L^4) )^0.5 ;
   oblique=exp( (betax1^2*alphax2+6*betax1*betax2/L^2+betax2^2*alphax1)/( 4*alphax1*alphax2-36/(L^4)) );
 
   fun3=w1*w2*pulse1*pulse2*perturb*oblique*pi^2/D3;
end   
