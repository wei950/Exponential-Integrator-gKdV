% 2nd order explicit low-regularity method for KdV; based on composition
% and 1st order prediction
% dispersionless limit: u_t+ep^2*u_xxx+1/k*(u^k)_x=0
nn=0;
N=2^9*4^nn; L=10*pi; h=2*L/N; x=-L:h:L-h; 
%k=1E-4/1;  
k=1/16;
T=100;
kini=round(T/k);  Nmax=kini;%k=T/kini;
ep=1/1;%1/4^nn;
mu=[0:N/2-1,-N/2:-1]*pi/L; emuh=exp(1i*mu.^3*k/2*ep^2); emu=exp(1i*mu.^3*k*ep^2);
%% initial data
nd=1; np=4; kp=np+1;  

%--rough data
s=4; 
u0=rand(1,N); u0=u0./abs(mu).^s; u0(1)=0; u0=ifft(u0); u0=u0/max(abs(u0)); uini=real(u0);
%u0=1.1*sin(x+pi);%uini(1:1:end);%
%--soliton data
% sigma=1;c=2; ep0=1E-1/2; x0=-.5;
ep0=1E-1/1;
% u0=uu;%1.4*(84*sech(5*x).^2).^(1/5);%2*sin(x+pi);%1.0*(tanh((x+2)/ep0)-tanh((x-2)/ep0))/2;%sigma*(c*(np+1)*(np+2)/2*(sech(sqrt(c)/ep*np/2*(x-x0))).^2).^(1/np);%1*sech(2*x.^2);% %sech(x/sqrt(ep)).^2;%sech(x).^2;sigma*(sech(K*(x-5))).^(2/np);%
% %track(1)=u0(N/2);
% uex=sigma*(c*(np+1)*(np+2)/2*(sech(sqrt(c)/ep*np/2*(x-x0-c*T))).^2).^(1/np);
ave=sum(real(u0).^(kp-1))*h/(2*L)*k/2; 
maxvalue=zeros(1,Nmax+1);maxvalue(1)=max(abs(u0));
 mass=zeros(Nmax+1,1); energy_4=mass;
 mass(1)=sqrt(sum(abs(u0).^2)*h);
 energy_4(1)=sum(ep^2/2*abs(ifft(fft(u0).*mu*1i)).^2-u0.^(kp+1)/kp/(kp+1))*h;
%% Time iteration
for n=1:Nmax
    %--1st order prediction    
    wk=u0.^(kp-1); wtemp=fft(u0).*emu;
    wet=fft(wk)./(1i*mu); wet(1)=0; 
    wt=fft(ifft(wet).*u0).*emu./(1i*mu); wt(1)=0; wk=wtemp+nd/3/ep^2*wt;
    wt=fft(ifft(wet.*emu).*ifft(wtemp))./(1i*mu); wt(1)=0; u=ifft(wk-nd/3/ep^2*wt); %prediction 
    
    wtemp=fft(u0).*emuh; wk=u0.^(kp-1);
    wet=fft(wk)./(1i*mu); wet(1)=0; 
    wt=fft(ifft(wet).*u0).*emuh./(1i*mu); wt(1)=0; wk=wtemp+nd/3/ep^2*wt;
    wt=fft(ifft(wet.*emuh).*ifft(wtemp))./(1i*mu); wt(1)=0; rhs=wk-nd/3/ep^2*wt;  
    
    wk=u.^(kp-1); wtemp=fft(u)./emuh;
    wet=fft(wk)./(1i*mu); wet(1)=0; 
    wt=fft(ifft(wet).*u)./emuh./(1i*mu); wt(1)=0; wk=nd/3/ep^2*wt;
    wt=fft(ifft(wet./emuh).*ifft(wtemp))./(1i*mu); wt(1)=0; u0=ifft((rhs-wk+nd/3/ep^2*wt).*emuh); 
    
    ave=ave+sum(real(u0).^(kp-1))*h/(2*L)*k/1;
      maxvalue(n+1)=max(abs(u0));
%        if mod(n,500)==0
%             u=fft(u0)/N*exp(1i*mu.'*(x+L-nd*ave));%track(n+1)=u(N/2);
%              plot(x,real(u)); %axis([-1,3,0,3]); 
%              drawnow
%        end
     mass(n+1)=sqrt(sum(abs(u0).^2)*h);
     energy_4(n+1)=sum(ep^2/2*abs(ifft(fft(u0).*mu*1i)).^2-u0.^(kp+1)/kp/(kp+1))*h;
end
save('mass_4.mat','mass_4')
save('energy_4.mat','energy_4')

%ave=ave-sum(real(u0).^(kp-1))*h/(2*L)*k/2;
%u1=fft(u0)/N*exp(1i*mu.'*(x+L-nd*ave));
% plot(x,real(u));hold on
%u0=fft(u0)/N*exp(1i*mu.'*(x+L-nd*ave)); %uex=u0;
%sqrt(sum(abs(ifft(fft(u0).*mu*1i)).^2)*h)
%eu=sqrt( sum( abs(uex(1:2^5:end)-u0).^2 ) )/sqrt( sum( abs(uex).^2 ) )

%plot(x,real(u0));%axis([0,3,0,3]);
%plot(x,u0)

%sqrt(sum(abs(u0-uex(1:1:end)).^2)*h)/sqrt(sum(abs(uex(1:1:end)).^2)*h)
plot(0.03:k:0.05,maxvalue)
hold on
%axis([0,T,0,100])
