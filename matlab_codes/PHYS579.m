clc; clear variables
% Pauli Matrices
sgma1 = [0 1;1 0];
sgma2 = [0 -1i;1i 0];
sgma3 = [1 0; 0 -1];
sgma0 = [1 0;0 1];

% Gamma matrices;

gma13 = kron(sgma1,sgma3);
gma30 = kron(sgma3,sgma0);
gma03 = kron(sgma0,sgma3);
gma12 = kron(sgma1,sgma2);
gma31 = kron(sgma3,sgma1);
gma21 = kron(sgma2,sgma1);
gma32 = kron(sgma3,sgma2);

NN = 28;
% Permutation
v = 1:NN;
C = nchoosek(v,NN/2);
ss = size(C);
t=1;
for ii=1:ss(1)
    for i = 1:ss(2)
        if C(ii,i) ~= 2*i && C(ii,i) ~= 2*i-1
            tt(t) = ii;
            t =t+1;
        end
    end
end

C(tt,:) = [];

m11 = C(1:2:length(C)/2,:);
m12 = C(2:2:length(C)/2,:);
m21 = C(length(C)/2+1:2:length(C),:);
m22 = C(length(C)/2+2:2:length(C),:);

% Variables
alpha1 = -1.5; beta1 = 1.5; gamma1 = 1;
% alpha1 = -1.5; beta1 = 0; gamma1 = 0;
% alpha1 = -1.5; beta1 = 1.5; gamma1 = -1;
kky = linspace(-pi,pi,250);
kkx  = linspace(0,2*pi,NN/2); 
wlp11 = zeros(1,length(kky)); wlp12 = zeros(1,length(kky)); wlp21 = zeros(1,length(kky)); wlp22 = zeros(1,length(kky));
for ii = 1:length(kky)
    for i = 1:length(kkx)
        hm = gma13/2 + (alpha1/2).*(cos(kkx(i)) + cos(kky(ii))).*(gma30 + gma03) + (gma12 + gma31).*sin(kkx(i)) + (gma21 + gma32).*sin(kky(ii)) - gma03 + (beta1/2).*cos(kkx(i)).*(cos(kky(ii))-gamma1).*(gma03 - gma30);
%         hm1 = gma13/2 + (alpha1/2).*(cos(kkx(i-1)) + cos(kky(ii))).*(gma30 + gma03) + (gma12 + gma31).*sin(kkx(i-1)) + (gma21 + gma32).*sin(kky(ii)) - gma03 + (beta1/2).*cos(kkx(i-1)).*(cos(kky(ii))-gamma1).*(gma03 - gma30);
        % --
        [V,D]=eig(hm);
        [d,ind] = sort(diag(D));
        VVs(:,2*i-1) = V(:,ind(1));
        VVs(:,2*i) = V(:,ind(2));
%         Vs = V(:,ind);
%         [V1,D1]=eig(hm1);
%         [d,ind] = sort(diag(D1));
%         Ds1 = D1(ind,ind);
%         Vs1 = V1(:,ind);
%         wlp11(ii) = wlp(ii)*(dot(Vs(:,1),Vs1(:,1))+dot(Vs(:,1),Vs1(:,2)));
%         wlp22(ii) = wlp1(ii)*(dot(Vs(:,2),Vs1(:,2))+dot(Vs(:,2),Vs1(:,1)));
    end

    for jj =1:length(m11(:,1))
        smwl11 = 1; smwl12 = 1; smwl21 = 1; smwl22 = 1;
        for j =1:length(m11(1,:))-1
            smwl11 = smwl11*dot(VVs(:,m11(jj,j)),VVs(:,m11(jj,j+1)));
            smwl12 = smwl12*dot(VVs(:,m12(jj,j)),VVs(:,m12(jj,j+1)));
            smwl21 = smwl21*dot(VVs(:,m21(jj,j)),VVs(:,m21(jj,j+1)));
            smwl22 = smwl22*dot(VVs(:,m22(jj,j)),VVs(:,m22(jj,j+1)));
        end
        wlp11(ii) = wlp11(ii) + smwl11;
        wlp12(ii) = wlp12(ii) + smwl12;
        wlp21(ii) = wlp21(ii) + smwl21;
        wlp22(ii) = wlp22(ii) + smwl22;
    end
    WilsonM = [wlp11(ii) wlp12(ii);wlp21(ii) wlp22(ii)];
    [Vw,Dw]=eig(WilsonM);
    e1wlp(ii) = Dw(1,1);
    e2wlp(ii) = Dw(2,2);
end
wnr = log(e1wlp)/(1i*2*pi);
wnr2 = log(e2wlp)/(1i*2*pi);

% A = [kky; wnr]; A2 = [kky; wnr2]; AA = A'; AA2 = A2';
% save('Wannier1.txt','AA','-ascii','-double'); save('Wannier12.txt','AA2','-ascii','-double');




%% plots
figure; hold on;
plot(real(wnr2),kky/(pi),'.')
plot(real(wnr),kky/(pi),'.')
plot(real(wnr2)+2,kky/(pi),'.')
plot(real(wnr)+2,kky/(pi),'.')
plot(real(wnr2)-2,kky/(pi),'.')
plot(real(wnr)-2,kky/(pi),'.')
plot(real(wnr2),(kky+2*pi)/pi,".")
plot(real(wnr),(kky+2*pi)/(pi),'.')
plot(real(wnr2)+2,(kky+2*pi)/(pi),'.')
plot(real(wnr)+2,(kky+2*pi)/(pi),'.')
plot(real(wnr2)-2,(kky+2*pi)/(pi),'.')
plot(real(wnr)-2,(kky+2*pi)/(pi),'.')
plot(real(wnr2),(kky-2*pi)/pi,".")
plot(real(wnr),(kky-2*pi)/(pi),'.')
plot(real(wnr2)+2,(kky-2*pi)/(pi),'.')
plot(real(wnr)+2,(kky-2*pi)/(pi),'.')
plot(real(wnr2)-2,(kky-2*pi)/(pi),'.')
plot(real(wnr)-2,(kky-2*pi)/(pi),'.')

%%
clear variables
v = 1:28;
C = nchoosek(v,14);
ss = size(C);
t=1;
for ii=1:ss(1)
    for i = 1:ss(2)
        if C(ii,i) ~= 2*i && C(ii,i) ~= 2*i-1
            tt(t) = ii;
            t =t+1;
        end
    end
end


%%


clc; clear variables
% Pauli Matrices
sgma1 = [0 1;1 0];
sgma2 = [0 -1i;1i 0];
sgma3 = [1 0; 0 -1];
sgma0 = [1 0;0 1];

% Gamma matrices;

gma13 = kron(sgma1,sgma3);
gma30 = kron(sgma3,sgma0);
gma03 = kron(sgma0,sgma3);
gma12 = kron(sgma1,sgma2);
gma31 = kron(sgma3,sgma1);
gma21 = kron(sgma2,sgma1);
gma32 = kron(sgma3,sgma2);

NN = 1000;

% Variables
alpha1 = -1.5; beta1 = 1.5; gamma1 = 1;
% alpha1 = -1.5; beta1 = 1.5; gamma1 = -1; NN = 1000; 
% alpha1 = -1.5; beta1 = 0; gamma1 = 0;
kky = pi/2; %linspace(-pi,pi,250);
kkx  = linspace(0,2*pi,NN); 
wlp11 = zeros(1,length(kky)); wlp21 = zeros(1,length(kky)); 
for ii = 1:length(kky)
    for i = 2:length(kkx)-1
        hm = gma13/2 + (alpha1/2).*(cos(kkx(i)) + cos(kky(ii))).*(gma30 + gma03) + (gma12 + gma31).*sin(kkx(i)) + (gma21 + gma32).*sin(kky(ii)) - gma03 + (beta1/2).*cos(kkx(i)).*(cos(kky(ii))-gamma1).*(gma03 - gma30);
        hm1 = gma13/2 + (alpha1/2).*(cos(kkx(i-1)) + cos(kky(ii))).*(gma30 + gma03) + (gma12 + gma31).*sin(kkx(i-1)) + (gma21 + gma32).*sin(kky(ii)) - gma03 + (beta1/2).*cos(kkx(i-1)).*(cos(kky(ii))-gamma1).*(gma03 - gma30);
        % --
        [V,D]=eig(hm);
        [d,ind] = sort(diag(D));
        Vs = V(:,ind);
        [V1,D1]=eig(hm1);
        [d,ind] = sort(diag(D1));
        Ds1 = D1(ind,ind);
        Vs1 = V1(:,ind);
        wlp11(i-1) = dot(Vs(:,1),Vs1(:,1));
        wlp22(i-1) = dot(Vs(:,2),Vs1(:,2));
    end
end

%% Assignment 2
clear variables
t0=zeros(9,9);
T0=zeros(6,6);
for i=1:length(T0)-1
    T0(i,i+1)=1;
end
for i=1:length(t0)
    t0(i,i)=1;
end
T0(length(T0),1)=1;
h0=kron(t0,T0+T0');

t=linspace(0,1.5,7);
Th1=zeros(6,6);
Th2=zeros(6,6);
Th3=zeros(6,6);
th1=zeros(9,9);
th2=zeros(9,9);
th3=zeros(9,9);
th1(1,2)=1; th1(2,4)=1; th1(3,5)=1; th1(5,7)=1; th1(6,8)=1; th1(8,9)=1;
th2(1,3)=1; th2(2,5)=1; th2(3,6)=1; th2(4,7)=1; th2(5,8)=1; th2(7,9)=1;
th3(2,3)=1; th3(4,5)=1; th3(5,6)=1; th3(7,8)=1;
for i=1:length(t)
    Th1=zeros(6,6);
    Th2=zeros(6,6);
    Th3=zeros(6,6);
    Th1(5,2)=t(i); Th2(4,1)=2*t(i); Th3(3,6)=t(i)/3;
    hh=kron(th1,Th1)+kron(th2,Th2)+kron(th3,Th3)+ctranspose(kron(th1,Th1))+ctranspose(kron(th2,Th2))+ctranspose(kron(th3,Th3));
    [V,D]=eig(h0+hh);
    En(i,:)=diag(D);
end
figure; hold on;
for i=1:54
    plot(t,En(:,i))
end

r(1,:)=[-0.6930603024962, 9.7009591155821];
r(2,:)=[0.6445441535552, 9.6648076437969];
r(3,:)=[1.439876532829, 8.5441120184566];
r(4,:)=[0.60839268177, 7.3149619777607];
r(5,:)=[-0.7292117742814, 7.2788105059755];
r(6,:)=[-1.5245441535552, 8.435657603101];


r(7,:)=[-3.1875118556731, 5.3627825013614];
r(8,:)=[-1.7776044560514, 5.4350854449317];
r(9,:)=[-1.0907264921331, 4.242086876021];
r(10,:)=[-1.7414529842662, 3.0129368353252];
r(11,:)=[-3.1513603838879, 2.97678536354];
r(12,:)=[-3.8382383478062, 4.2782383478062];


r(13,:)=[1.7290883071103, 5.4350854449317];
r(14,:)=[3.1028442349469, 5.4350854449317];
r(15,:)=[3.8258736706503, 4.2782383478062];
r(16,:)=[3.1028442349469, 3.0490883071103];
r(17,:)=[1.7290883071103, 3.1575427224659];
r(18,:)=[1.0422103431921, 4.242086876021];


r(19,:)=[-5.5735089934945, 1.2053632460666];
r(20,:)=[-4.1997530656579, 1.2776661896369];
r(21,:)=[-3.4405721581693, 0.0123646771559];
r(22,:)=[-4.1997530656579, -1.21678536354];
r(23,:)=[-5.6458119370648, -1.21678536354];
r(24,:)=[-6.3688413727683, -0.0960897381997];


r(25,:)=[-0.7653632460666, 1.3138176614221];
r(26,:)=[0.6445441535552, 1.2415147178517];
r(27,:)=[1.4037250610438, 0.0123646771559];
r(28,:)=[0.6806956253403, -1.1806338917548];
r(29,:)=[-0.6930603024962, -1.1806338917548];
r(30,:)=[-1.4522412099848, 0.048516148941];


r(31,:)=[4.1150854449317, 1.2776661896369];
r(32,:)=[5.633447259909, 1.3138176614221];
r(33,:)=[6.3564766956124, 0.0123646771559];
r(34,:)=[5.5972957881238, -1.2529368353252];
r(35,:)=[4.1873883885021, -1.2529368353252];
r(36,:)=[3.536661896369, -0.0599382664145];


r(37,:)=[-3.2236633274583, -2.9159045374431];
r(38,:)=[-1.6691500406959, -2.9882074810135];
r(39,:)=[-1.0184235485628, -4.145054578139];
r(40,:)=[-1.7776044560514, -5.3742046188348];
r(41,:)=[-3.2236633274583, -5.3019016752645];
r(42,:)=[-3.8382383478062, -4.2535089934945];


r(43,:)=[1.8013912506807, -2.9159045374431];
r(44,:)=[3.1028442349469, -3.0243589527986];
r(45,:)=[3.8258736706503, -4.3258119370648];
r(46,:)=[3.1389957067321, -5.3742046188348];
r(47,:)=[1.8013912506807, -5.3742046188348];
r(48,:)=[1.1145132867624, -4.145054578139];


r(49,:)=[-0.656908830711, -7.3625355670193];
r(50,:)=[0.7168470971255, -7.3263840952341];
r(51,:)=[1.4037250610438, -8.4109282487893];
r(52,:)=[0.7529985689107, -9.6400782894852];
r(53,:)=[-0.6207573589259, -9.7123812330555];
r(54,:)=[-1.3076353228441, -8.4470797205745];

x=linspace(-8,8,100); y = linspace(-11,11,100);
[X,Y]=meshgrid(x,y);
Z=0;
for i=1:length(diag(D))
    Z=Z+10*(V(i,27)-V(i,28)).*exp(-((X-r(i,1)).^2./(0.1)+(Y-r(i,2)).^2./(0.1)));
%     Z=Z+10*(V(i,27)).*exp(-((X-r(i,1)).^2./(0.1)+(Y-r(i,2)).^2./(0.1)));
end
figure;
s = mesh(X,Y,Z,'FaceAlpha','0.5'); hold on
axis equal

%% 1d)

clear variables

Gamma = [0, 0];
K = [2*pi/3 0];
KK = [pi/3 sqrt(3)*pi/3];
M = [pi/2 pi/(2*sqrt(3))];
t = 1.5;
alpha1 = [2 0];
alpha2 = [1 sqrt(3)];
alpha3 = [1 -sqrt(3)];
h1 = [0 1 0 t*exp(-1i*dot(Gamma,alpha3)) 0 1;1 0 1 0 t*exp(1i*dot(Gamma,alpha2)) 0;0 1 0 1 0 t*exp(1i*dot(Gamma,alpha1));t*exp(1i*dot(Gamma,alpha3)) 0 1 0 1 0;0 t*exp(-i*dot(Gamma,alpha2)) 0 1 0 1;1 0 t*exp(-1i*dot(Gamma,alpha1)) 0 1 0];
h2 = [0 1 0 t*exp(-1i*dot(K,alpha3)) 0 1;1 0 1 0 t*exp(1i*dot(K,alpha2)) 0;0 1 0 1 0 t*exp(1i*dot(K,alpha1));t*exp(1i*dot(K,alpha3)) 0 1 0 1 0;0 t*exp(-i*dot(K,alpha2)) 0 1 0 1;1 0 t*exp(-1i*dot(K,alpha1)) 0 1 0];
h3 = [0 1 0 t*exp(-1i*dot(KK,alpha3)) 0 1;1 0 1 0 t*exp(1i*dot(KK,alpha2)) 0;0 1 0 1 0 t*exp(1i*dot(KK,alpha1));t*exp(1i*dot(KK,alpha3)) 0 1 0 1 0;0 t*exp(-i*dot(KK,alpha2)) 0 1 0 1;1 0 t*exp(-1i*dot(KK,alpha1)) 0 1 0];
h4 = [0 1 0 t*exp(-1i*dot(M,alpha3)) 0 1;1 0 1 0 t*exp(1i*dot(M,alpha2)) 0;0 1 0 1 0 t*exp(1i*dot(M,alpha1));t*exp(1i*dot(M,alpha3)) 0 1 0 1 0;0 t*exp(-i*dot(M,alpha2)) 0 1 0 1;1 0 t*exp(-1i*dot(M,alpha1)) 0 1 0];
R2=zeros(6,6);
R2(1,4)=1; R2(2,5)=1; R2(3,6)=1; R2(4,1)=1; R2(5,2)=1; R2(6,3)=1;
R3=zeros(6,6);
R3(1,3)=1; R3(2,4)=1; R3(3,5)=1; R3(4,6)=1; R3(5,1)=1; R3(6,2)=1;
[V,D]=eig(h2);

for j=1:2
    for jj=1:2
        hmd(j,jj)=(V(:,j)')*R3*V(:,jj);
    end
end