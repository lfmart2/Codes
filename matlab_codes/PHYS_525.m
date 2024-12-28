clc; clear all;

t=1; a=1; kx=linspace(-pi/a,pi/a,1000); ky=linspace(-pi/a,pi/a,1000);
DelE=linspace(0,8,8000);Energy=-4; error=0.05;
for i=1:length(DelE)
    rho(i)=0;
    for ii=1:length(kx)
        for j=1:length(ky)
            prf=Energy+DelE(i)-2*t*(cos(kx(ii)*a)+cos(ky(j)*a));
            if abs(prf) < error
                rho(i)=rho(i)+1;
            end
        end
    end
    rho(i)=rho(i)/(2*pi^2);
end
figure; hold on;
plot(Energy+DelE,rho)
ylabel('DOS $\rho(E)$','Interpreter','latex')
xlabel('Energy $E(\vec{k})$','Interpreter','latex')

%%
clc; clear all;

t=1; a=1; k=linspace(-pi/a,pi/a,1000);
DelE=linspace(0,8,8000);Energy=-1; error=0.5; v=10; w=10;
for i=1:length(DelE)
    rho(i)=0;
    for j=1:length(k)
        prf=Energy+DelE(i)-sqrt(v^2+w^2+2*w*v*cos(k(j)));
        if abs(prf) < error
            rho(i)=rho(i)+1;
        end
    end
    rho(i)=rho(i)/(2*pi^2);
end
figure; hold on;
plot(Energy+DelE,rho)
ylabel('DOS $\rho(E)$','Interpreter','latex')
xlabel('Energy $E(\vec{k})$','Interpreter','latex')