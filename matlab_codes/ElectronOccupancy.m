%% Electron density in the EdgeModes
clc; clear all;
v=10; w=linspace(5,15,400);
e0=linspace(0,0.1,5); 
N=800; T1=0.1;
n0 = zeros([length(e0) length(w)]);
n01 = zeros([length(e0) length(w)]);
for i=2:1:length(e0)
    for jj=1:length(w)
        psi=SSH_chain(N,v,w(jj));
        psi(1,2*N+1)=T1; psi(2*N+1,1)=T1;
        psi(2,2*N+1)=T1/3; psi(2*N+1,2)=T1/3;
        psi(2*N+1,2*N+1) = e0(i);
        [V,D]=eig(psi);
        Enner(jj,:)=diag(D);
        n0(i,jj) = sum(V(2*N+1, 1:N).^2, 'all');
        psi=SSH_chain(N,v,w(jj));
        [V,D]=eig(psi);
        n01(i,jj) = ((T1*max( abs( V(1,:) ) )/e0(i) )^2)/(1+(T1*max( abs( V(1,:) ) )/e0(i) )^2);
        n02(i,jj) = ((T1*max( abs( V(1,:) ) )/e0(i) )^2);
    end
end
figure; hold on;   
plot(w./v,n0(1,:),'LineWidth',2);
plot(w./v,n0(2,:),'--','LineWidth',2);
plot(w./v,n0(3,:),':','LineWidth',2);
plot(w./v,n0(4,:),'-.','LineWidth',2);
plot(w./v,n0(5,:),'LineWidth',2);
plot(w./v,n01(1,:),'LineWidth',2);
plot(w./v,n01(2,:),'--','LineWidth',2);
plot(w./v,n01(3,:),':','LineWidth',2);
plot(w./v,n01(4,:),'-.','LineWidth',2);
plot(w./v,n01(5,:),'LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlim([w(1)/v w(end)/v])
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
box on;
l = legend('show','interpreter','latex','FontWeight','bold');
l.Title.String = '$\bf{\varepsilon_0}$';
l.LineWidth = 0.5;
fontSize = 15;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
legend('$0$','$0.025$','$0.05$','$0.075$','0.1','interpreter','latex');
ylim([0 0.65])

%% Electron density in the bulk
clear all; clc
tic
v=10; w=linspace(10.0,10.06,500); %10.375;
T1=0.1;[0.05 0.1 0.15]; N=2*800; e0=0;%[0 4*10^-4 8*10^-4 12*10^-4];
n0 = zeros(length(T1),length(w));
for i=1:1:length(T1)
    for jj=1:length(w)
        psi=SSH_chain(N,v,w(jj));
        % psi=SSH(N,v,w(jj),e0(i));
        psi(2*N+1,2*N+1) = e0;
        psi(N,2*N+1)=T1(i); psi(2*N+1,N)=T1(i); 
        psi(N+1,2*N+1)=T1(i)/3; psi(2*N+1,N+1)=T1(i)/3;
        psi(N-1,2*N+1)=T1(i)/3; psi(2*N+1,N-1)=T1(i)/3;
        [V,D]=eig(psi);
        Enner(jj,:)=diag(D);
        % 
        % [~, maxIndex] = max(abs(V(2*N+1,:)));
        % V(:,maxIndex) = [];
        % n0(i,jj) = sum(V(2*N+1, 1:N).^2, 'all');
        n0(i,jj) = max(abs(V(1,:)));
         
        % if abs(V(2*N+1,N)) < 0.65
        %     for ii=1:N 
        %         n0(i,jj)=n0(i,jj)+(V(length(psi),ii))^2;
        %     end
        % else
        %     for ii=1:N-1
        %         n0(i,jj)=n0(i,jj)+(V(length(psi),ii))^2;
        %     end
        %     n0(i,jj)=n0(i,jj)+(V(length(psi),N+1))^2;
        % end
    end
end
figure; hold on;   
plot(w./v,n0(1,:),'LineWidth',2);
% plot(w./v,n0(2,:),'--','LineWidth',2);
% plot(w./v,n0(3,:),':','LineWidth',2);
% plot(w./v,n0(4,:),'-.','LineWidth',2);
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_B)$','interpreter','latex')
l = legend('show','interpreter','latex','FontWeight','bold');
% l.Title.String = '$\bf{\varepsilon_0}$';
l.Title.String = '$\bf{T_{m,\alpha}}$';
legend('$0.05$','$0.1$','$0.15$','interpreter','latex');
% legend('$0$','$0.0004$','$0.0008$','$0.0012$','interpreter','latex');
xlim([w(1)/v w(end)/v])
fontSize = 15;
box on;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
toc



%% Comparison of edge Modes and bulk
clc; clear all;
v=10; w=15;
e0=linspace(0,0.15,400); N=800; T1=linspace(0.05,0.15,3);
figure; hold on
for i=1:1:length(T1)
    for jj=1:length(e0)
        psi=SSH_chain(N,v,w);
        psi(1,2*N+1)=T1(i); psi(2*N+1,1)=T1(i);
        psi(2,2*N+1)=T1(i)/3; psi(2*N+1,2)=T1(i)/3;
        psi(2*N+1,2*N+1) = e0(jj);
        [V,D]=eig(psi);
        n0(i,jj)=0;       
        for ii=1:N %length(psi)
            n0(i,jj)=n0(i,jj)+(V(length(psi),ii))^2;
        end
    end
end
vv=10; ww=[10.0814 10.0597 10.0465];
for i=1:1:length(T1)
    for jj=1:1:length(e0)
        psi=SSH_chain(N,vv,ww(i));
        psi(N,2*N+1)=-T1(i); psi(2*N+1,N)=-T1(i); 
        psi(N+1,2*N+1)=-T1(i)/3; psi(2*N+1,N+1)=-T1(i)/3;
        psi(N-1,2*N+1)=-T1(i)/3; psi(2*N+1,N-1)=-T1(i)/3;
        psi(2*N+1,2*N+1) = e0(jj);
        [V1,D]=eig(psi);
        n00(i,jj)=0;
        % if abs(V1(2*N+1,N)) < 0.60
        %     for ii=1:N %length(psi)
        %         n00(i,jj)=n00(i,jj)+(V1(length(psi),ii))^2;
        %     end
        % else
        %     for ii=1:N-1 %length(psi)
        %         n00(i,jj)=n00(i,jj)+(V1(length(psi),ii))^2;
        %     end
        %     n00(i,jj)=n00(i,jj)+(V1(length(psi),N+1))^2;
        % end
        [~, maxIndex] = max(abs(V1(2*N+1,:)));
        V1(:,maxIndex) = [];
        n00(i,jj) = sum(V1(2*N+1, 1:N).^2, 'all');
    end
end
xlabel('$\varepsilon_0$','interpreter','latex'); ylabel('$n_0(R,x_M)$','interpreter','latex')
% title('$N=200$', 'interpreter','latex');
l = legend('show','interpreter','latex','FontWeight','bold');
l.Title.String = '$T_{1,A},\, T_{m,A}$';
% title('$N=200$, $v=10$, and $T_1=3T_3=0.1$','interpreter','latex');
xlim([0 0.15]); ylim([0 0.5]);
fontSize = 15;
box on;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);

blueColor = [0, 0.4470, 0.7410];
redColor = [0.85,0.33,0.10];
yellowColor = [0.93,0.69,0.13];

plot(e0,n0(1,:),'Color',blueColor,'LineWidth',2,'HandleVisibility','off');
plot(e0,n0(2,:),'Color',redColor,'LineWidth',2,'HandleVisibility','off');
plot(e0,n0(3,:),'Color',yellowColor,'LineWidth',2,'HandleVisibility','off'); 

plot(e0,n00(1,:),'--','Color',blueColor,'LineWidth',2,'HandleVisibility','off');
plot(e0,n00(2,:),'--','Color',redColor,'LineWidth',2,'HandleVisibility','off');
plot(e0,n00(3,:),'--','Color',yellowColor,'LineWidth',2,'HandleVisibility','off'); 

x = [0.1]; y = [0.6]; 
plot(x,y,'o','Color',blueColor,'MarkerFaceColor',blueColor)
plot(x,y,'o','Color',redColor,'MarkerFaceColor',redColor)
plot(x,y,'o','Color',yellowColor,'MarkerFaceColor',yellowColor)
legend('$0.05$','$0.1$','$0.15$','interpreter','latex');
%% Bulk Defect
clc; clear all;
v=10; w=linspace(5,15,400);
T1 = 0.1; e0=[0.025 0.075 0.1]; 
N=800;
for j=1:length(e0)
    for jj=1:length(w)        
        [V,D]=eig(DSSH(N,T1,v,w(jj),e0(j),e0(j)));
        Enner(jj,:)=diag(D);
        n0(j,jj)=0;

        n00(j,jj)=0;
        for ii=1:N/2+1
            n0(j,jj)=n0(j,jj)+(V(N+2,ii))^2;
            n00(j,jj)=n00(j,jj)+(V(N+3,ii))^2;
        end
%         n0(j,jj)=n0(j,jj)-(V(2*N+2,N+2))^2;
        ntl(j,jj) = n0(j,jj) + n00(j,jj);
    end
end


% plot(w./v,n0(1,:)); 
% plot(w./v,n0(2,:)); plot(w./v,n0(3,:))
% xline(1,'-.')
% xlabel('$w/v$','Interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
% legend('$\epsilon_0(R_1)=0.025$','$\epsilon_0(R_3)=0.075$','$\epsilon_0(R_4)=0.1$','interpreter','latex');
% title('Interaction in the edges','interpreter','latex');
% 
% figure;  hold on;
% plot(w./v,n00(1,:));
% plot(w./v,n00(2,:)); plot(w./v,n00(3,:));
% xline(1,'-.')
% xlabel('$w/v$','Interpreter','latex'); ylabel('$n_0(R,x_B)$','interpreter','latex')
% legend('$\epsilon_0(R_1)=0.025$','$\epsilon_0(R_3)=0.075$','$\epsilon_0(R_4)=0.1$','interpreter','latex');
% title('Interaction in the bulk','interpreter','latex');
% figure;  hold on;
% 
% plot(w./v,ntl(1,:)); plot(w./v,ntl(2,:)); plot(w./v,ntl(3,:));
% xline(1,'-.')
% xlabel('$w/v$','Interpreter','latex'); ylabel('$\Delta n_0(R,x_B)$','interpreter','latex')
% legend('$\epsilon_0(R_1)=0.025$','$\epsilon_0(R_3)=0.075$','$\epsilon_0(R_4)=0.1$','interpreter','latex');
% title('Total occupation number','interpreter','latex');

figure; hold on
subplot(3,3,[1,2,3])
plot(w./v,n0(1,:),'LineWidth',2); hold on;
plot(w./v,n0(2,:),'--','LineWidth',2); plot(w./v,n0(3,:),':','LineWidth',2)
xline(1,'-.','HandleVisibility','off')
ylabel('$n_0(R,x_E)$','interpreter','latex')
title('Interaction in the edges','interpreter','latex');
l = legend('show','interpreter','latex','FontWeight','bold');
l.Title.String = '$\bf{\varepsilon_0}$';
legend('$0.025$','$0.075$','$0.1$','interpreter','latex');
fontSize = 15;
box on;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);

subplot(3,3,[4,5,6])
plot(w./v,n00(1,:),'LineWidth',2); hold on;
plot(w./v,n00(2,:),'--','LineWidth',2); plot(w./v,n00(3,:),':','LineWidth',2);
xline(1,'-.','HandleVisibility','off')
ylabel('$n_0(R,x_B)$','interpreter','latex')
title('Interaction in the bulk','interpreter','latex');
box on;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);

subplot(3,3,[7,8,9])
plot(w./v,ntl(1,:),'LineWidth',2); hold on
plot(w./v,ntl(2,:),'--','LineWidth',2); plot(w./v,ntl(3,:),':','LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlabel('$w/v$','Interpreter','latex'); ylabel('$\Sigma n_0(R,x_B)$','interpreter','latex')
title('Total occupation number','interpreter','latex');
box on;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);

%% Plots for comparison:
figure; hold on;
subplot(2,2,1)
load('met_N800.mat')
plot(w./v,n0(1,:),'LineWidth',2);
hold on;
title('N = 800')
plot(w./v,n0(2,:),'--','LineWidth',2);
plot(w./v,n0(3,:),':','LineWidth',2);
plot(w./v,n0(4,:),'-.','LineWidth',2);
plot(w./v,n0(5,:),'LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlim([w(1)/v w(end)/v])
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
% box on;
% l = legend('show','interpreter','latex','FontWeight','bold');
% l.Title.String = '$\bf{\varepsilon_0}$';
% l.LineWidth = 0.5;
fontSize = 15;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
% legend('$0$','$0.025$','$0.05$','$0.075$','0.1','interpreter','latex');
ylim([0 0.65])

subplot(2,2,2)
load('met_N1000.mat')
plot(w./v,n0(1,:),'LineWidth',2);
hold on;
title('N = 1000')
plot(w./v,n0(2,:),'--','LineWidth',2);
plot(w./v,n0(3,:),':','LineWidth',2);
plot(w./v,n0(4,:),'-.','LineWidth',2);
plot(w./v,n0(5,:),'LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlim([w(1)/v w(end)/v])
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
box on;
l = legend('show','interpreter','latex','FontWeight','bold');
l.Title.String = '$\bf{\varepsilon_0}$';
l.LineWidth = 0.5;
fontSize = 15;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
e0=[0 10^-4 10^-3 0.05 0.1];
legend('$0$','$10^{-4}$','$10^{-3}$','$0.05$','0.1','interpreter','latex');
ylim([0 0.65])


subplot(2,2,3)
load('met_N1500.mat')
plot(w./v,n0(1,:),'LineWidth',2);
hold on; title('N = 1500')
plot(w./v,n0(2,:),'--','LineWidth',2);
plot(w./v,n0(3,:),':','LineWidth',2);
plot(w./v,n0(4,:),'-.','LineWidth',2);
plot(w./v,n0(5,:),'LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlim([w(1)/v w(end)/v])
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
% box on;
% l = legend('show','interpreter','latex','FontWeight','bold');
% l.Title.String = '$\bf{\varepsilon_0}$';
% l.LineWidth = 0.5;
fontSize = 15;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
% legend('$0$','$0.025$','$0.05$','$0.075$','0.1','interpreter','latex');
ylim([0 0.65])


subplot(2,2,4)
load('met_N2000.mat')
plot(w./v,n0(1,:),'LineWidth',2);
hold on;
title('N = 2000')
plot(w./v,n0(2,:),'--','LineWidth',2);
plot(w./v,n0(3,:),':','LineWidth',2);
plot(w./v,n0(4,:),'-.','LineWidth',2);
plot(w./v,n0(5,:),'LineWidth',2);
xline(1,'-.','HandleVisibility','off')
xlim([w(1)/v w(end)/v])
xlabel('$w/v$','interpreter','latex'); ylabel('$n_0(R,x_E)$','interpreter','latex')
% box on;
% l = legend('show','interpreter','latex','FontWeight','bold');
% l.Title.String = '$\bf{\varepsilon_0}$';
% l.LineWidth = 0.5;
fontSize = 15;
set(gca, 'FontSize', fontSize, 'LineWidth', 1);
% legend('$0$','$0.025$','$0.05$','$0.075$','0.1','interpreter','latex');
ylim([0 0.65])

%% Max points
clc; clear all
blueColor = [0, 0.4470, 0.7410];
redColor = [0.85,0.33,0.10];
yellowColor = [0.93,0.69,0.13];
purpleColor = [0.49,0.18,0.56];
greenColor = [0.4660 0.6740 0.1880];
load('met_N800.mat')
[val(1) idx(1)] = max(n0(1,:));
[val(2) idx(2)] = max(n0(2,:));
[val(3) idx(3)] = max(n0(3,:));
ww = w;
load('met_N1000.mat')
[val1(1) idx1(1)] = max(n0(1,:));
[val1(2) idx1(2)] = max(n0(2,:));
[val1(3) idx1(3)] = max(n0(3,:));
ww1 = w;
load('met_N1500.mat')
[val2(1) idx2(1)] = max(n0(1,:));
[val2(2) idx2(2)] = max(n0(2,:));
[val2(3) idx2(3)] = max(n0(3,:));
ww2 = w;
load('met_N2000.mat')
[val3(1) idx3(1)] = max(n0(1,:));
[val3(2) idx3(2)] = max(n0(2,:));
[val3(3) idx3(3)] = max(n0(3,:));
ww3 = w;
figure; hold on;
scatter(ww(idx(1))/v,val(1),'*','MarkerEdgeColor',blueColor);
scatter(ww(idx(2))/v,val(2),'+','MarkerEdgeColor',blueColor,'HandleVisibility','off');
scatter(ww(idx(3))/v,val(3),'x','MarkerEdgeColor',blueColor,'HandleVisibility','off');
scatter(ww1(idx1(1))/v,val1(1),'*','MarkerEdgeColor',yellowColor);
scatter(ww1(idx1(2))/v,val1(2),'+','MarkerEdgeColor',yellowColor,'HandleVisibility','off');
scatter(ww1(idx1(3))/v,val1(3),'x','MarkerEdgeColor',yellowColor,'HandleVisibility','off');

scatter(ww2(idx2(1))/v,val2(1),'*','MarkerEdgeColor',redColor);
scatter(ww2(idx2(2))/v,val2(2),'+','MarkerEdgeColor',redColor,'HandleVisibility','off');
scatter(ww2(idx2(3))/v,val2(3),'x','MarkerEdgeColor',redColor,'HandleVisibility','off');

scatter(ww3(idx3(1))/v,val3(1),'*','MarkerEdgeColor',purpleColor);
scatter(ww3(idx3(2))/v,val3(2),'+','MarkerEdgeColor',purpleColor,'HandleVisibility','off');
scatter(ww3(idx3(3))/v,val3(3),'x','MarkerEdgeColor',purpleColor,'HandleVisibility','off');
xlabel("$w/v$",'interpreter','latex')
ylabel('Electron occupancy')
legend('$N = 800$','$N = 1000$','$N = 1500$','$N = 2000$','interpreter','latex')