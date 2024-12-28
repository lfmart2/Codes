%% Soliton length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% calculation of the Soliton length l for a given v/w value %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
tic
% SSH Chain parameters
N= 10000; v = 0.59; w = 0.88; 

% SSH chain
[Vs DDs] = eig(SSH_chain(N/2,v,w));
d = diag(DDs);
d(length(d)/2+1:end) = [];
smd = sum(d);

long = linspace(100,150,50);

for i=1:length(long)
    ss1 = ((w-v)*tanh((linspace(1,N-1,N-1)-N/2)./long(i)))./2;
    ss1(1:2:length(ss1)) = -ss1(1:2:length(ss1));
    ss = (w+v)/2+ss1;
    
    solSSH = diag(ss,1) + ctranspose(diag(ss,1));
    
    dlt = zeros([1 length(ss)]);
    for ii = 1:length(ss)
        dlt(ii) = ss(ii)^2;
    end
    Ds = sort(eig(solSSH));
    Ds(length(Ds)/2+1:end) = [];
    En_sol(i) = sum(Ds) - smd + sum(dlt)-(N/2-1)*w - v*N/2;
end

toc
plot(long,En_sol)

%%
clear variables
tic
% SSH Chain parameters
N= 6078; v = 10; w = 15;%linspace(7,13,5);
T1 = 0.1; e0 = 0.1;
for jj=1:length(w)
    % % SSH chain
    % d = sort(eig(SSH_chain(N/2,v,w(jj))));
    % d(length(d)/2+1:end) = [];
    % smd = sum(d);
    long = 7;%linspace(0.01,10,50);
    
    for i=1:length(long)
        ss1 = (w(jj)-v)*(tanh((linspace(1,N-1,N-1)-N/3)./long(i))-tanh((linspace(1,N-1,N-1)-2*N/3)./long(i)))./2;
        ss1(2:2:length(ss1)) = -ss1(2:2:length(ss1));
        tp = (w(jj)-v)*ones(1,N-1)/2; tp(1:2:length(tp)) = -tp(1:2:length(tp));
        ss = (w(jj)+v)/2+tp + ss1;
        
        solSSH = diag(ss,1) + ctranspose(diag(ss,1));
        
        % Ds = sort(eig(solSSH));
        % Ds(length(Ds)/2+1:end) = [];
        % En_sol(jj,i) = sum(Ds) - smd;
    end
end
solSSH(1,N+1) = -T1; solSSH(N+1,1) = -T1;
solSSH(2,N+1) = -T1/3; solSSH(N+1,2) = -T1/3;
solSSH(N+1,N+1) = e0;
solSSH(N/3,N+1+1) = -T1; solSSH(N+1+1,N/3) = -T1;
solSSH(N/3+1,N+1+1) = -T1; solSSH(N+1+1,N/3+1) = -T1;
% solSSH(nn(i)-1,N+i+1) = -T1/3; solSSH(N+i+1,nn(i)-1) = -T1/3;
solSSH(N+1+1,N+1+1) = e0;
solSSH(2*N/3,N+2+1) = -T1; solSSH(N+2+1,2*N/3) = -T1;
solSSH(2*N/3+1,N+2+1) = -T1; solSSH(N+2+1,2*N/3+1) = -T1;
% solSSH(nn(i)-1,N+i+1) = -T1/3; solSSH(N+i+1,nn(i)-1) = -T1/3;
solSSH(N+2+1,N+2+1) = e0;

[V D] = eig(solSSH);
n00 = zeros([1 3]);
for n= 1:3
    for jj = 1:N/2
        n00(n) = n00(n) + (V(N+n,jj))^2;
    end
end
toc
% figure; hold on
% for i = 1:length(w)
%     plot(long,En_sol(i,:));
% end