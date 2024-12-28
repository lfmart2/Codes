%% Poster SETCA plot
subplot(4,4,[1,2,5,6]);
k=linspace(-pi,pi,100);
v = 15; w = 10;
plot(k,sqrt(v^2+w^2+2*w*v*cos(k))); hold on
plot(k,-sqrt(v^2+w^2+2*w*v*cos(k)))
xline(0,'-.')
yline(0,'-.')
ylabel('Energy','Interpreter','latex'); xlabel('$k$','Interpreter','latex')
title('Trivial Insulator phase: $v\gg w$','Interpreter','latex')

subplot(4,4,[3,4,7,8]);
k=linspace(-pi,pi,100);
v = 10; w = 10;
plot(k,sqrt(v^2+w^2+2*w*v*cos(k))); hold on
plot(k,-sqrt(v^2+w^2+2*w*v*cos(k)))
xline(0,'-.')
yline(0,'-.')
ylabel('Energy','Interpreter','latex'); xlabel('$k$','Interpreter','latex')
title('Metalic phase: v=w','Interpreter','latex')

subplot(4,4,[9,10,13,14]);
k=linspace(-pi,pi,100);
v = 10; w = 15;
plot(k,sqrt(v^2+w^2+2*w*v*cos(k))); hold on
plot(k,-sqrt(v^2+w^2+2*w*v*cos(k)))
xline(0,'-.')
yline(0,'-.')
ylabel('Energy','Interpreter','latex'); xlabel('$k$','Interpreter','latex')
title('Topological phase: $w\gg v$','Interpreter','latex')

subplot(4,4,[11,12]);
v = 10; w = 15; N=20;
psi=zeros(2*N);
psi(1,2)=v; psi(2*N,2*N-1)=v;
for ii=2:2:2*(N-1)
    psi(ii,ii-1)=v;
    psi(ii,ii+1)=w;
    psi(ii+1,ii)=w;
    psi(ii+1,ii+2)=v;
end
[V,D]=eig(psi);
% for i=1:2:length(V(:,N))
%     VV(i)=V(i,N);
%     V(i,N)=0;
%     VV1(i+1)=V(i+1,N);
%     V(i,N+1)=0;
% end
bar(-(V(:,N)+V(:,N+1))/sqrt(2));
% hold on
% bar(VV);
subplot(4,4,[15,16])
bar((V(:,N)-V(:,N+1))/sqrt(2)); hold on
bar((V(:,N)-V(:,N+1))/sqrt(2));
% hold on
% bar(VV1);


%% Metalic transition dependence on T1
clc; clear all;
v=10; w=linspace(9,11,401);
T1=linspace(0,1,4); N=200; e0=0;
figure
for i=1:1:length(T1)
    for jj=1:length(w)
        psi=zeros(2*N+1);
        psi(1,2)=v; psi(2*N,2*N-1)=v;
        psi(N,2*N+1)=-T1(i); psi(2*N+1,N)=-T1(i); 
%         psi(N+1,2*N+1)=-T1/3; psi(2*N+1,N+1)=-T1/3;
%         psi(N-1,2*N+1)=-T1/3; psi(2*N+1,N-1)=-T1/3;
        psi(2*N+1,2*N+1)=e0;
        for ii=2:2:2*(N-1)
            psi(ii,ii-1)=v;
            psi(ii,ii+1)=w(jj);
            psi(ii+1,ii)=w(jj);
            psi(ii+1,ii+2)=v;
        end
        [V,D]=eig(psi);
        Enner(jj,:)=diag(D);
    end
    subplot(2,2,i)
    plot(w./v,Enner(:,N))
    hold on
    plot(w./v,Enner(:,N+1))
    plot(w./v,Enner(:,N+2))
    xline(1)
    xlabel('$w/v$','interpreter','latex'); ylabel('$n_0$','interpreter','latex');
    title(['N=200, v=10, ε_0=0 and T_1=' num2str(T1(i))]);
end


%%
figure;
t = tiledlayout(1,2,'TileSpacing','compact');
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];
ax1 = axes(t);
bar(ax1,V(:,N+1))
xline(ax1,N+20,':');
ax1.Box = 'off';
xlim(ax1,[N-20 N+20])
xlabel(ax1,'Sublattice')
ax2 = axes(t);
ax2.Layout.Tile = 2;
bar(ax2,V(:,N+1))
ax2.YAxis.Visible = 'off';
xline(ax2,2*N-15,':');
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[2*N-15 2*N+3])
xlabel(ax2,'Sublattice')
linkaxes([ax1 ax2], 'y')
title(ax2,'Adparticle contribution','Interpreter','latex')
title(ax1,'Bulk Defect Mode','Interpreter','latex')
ylabel(ax1,'Electron Occupancy $n_0$ at $w/v = 2$','Interpreter','latex')