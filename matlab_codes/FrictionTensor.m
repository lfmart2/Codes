%% Friction Tensor in the edges
clear variables
tic
beta = 1/0.015; e0 = 0.15; v = 10; w = 15; g = 0.08;
ef = 0; N=800; T1A=0.15; T1B = T1A/3; psi=zeros(2*N);
edi = e0-0.45; edf = e0+0.15;
% --
%Construction of the unperturbed Hamiltonian:
psi(1,2)=v; psi(2*N,2*N-1)=v;
for ii=2:2:2*(N-1)
    psi(ii,ii-1)=v;
    psi(ii,ii+1)=w;
    psi(ii+1,ii)=w;
    psi(ii+1,ii+2)=v;
end
% --
% Diagonalization of the matrix
[V,D]=eig(psi);
% --
% Sort of the eigenstates and energies
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
% --
% Linewidth of the Deltas function will be between the finite
% difference between the two states with closest energy difference: the
% first and second state
clear y Gamma Vk; nM=0.4*(Ds(2,2)-Ds(1,1)); l = 1; chck = 0;
% --
% Number of steps taken in consideration for the construction of the
% Gamma and Lambda Function. This term is the more sensitive and a small
% change gives unsense results


% --
% Time steps taken to evaluate the integration
ed = g*sqrt(2).*linspace((edi-e0)/(sqrt(2)*g),(edf-e0)/(sqrt(2)*g),(edf-edi)/(1*(Ds(2,2)-Ds(1,1)))) + e0;
eps = linspace(Ds(1,1),Ds(2*N,2*N),(Ds(2*N,2*N)-Ds(1,1))/(1*(Ds(2,2)-Ds(1,1)))); 
% --
% The term V(1,:) gives the projection of |1_A> into the d_k^\dagger. 
% The other term V(:,1) gives the eigenstate with energy E_1

Vk1 = Vs(1,:); Vk2 = Vs(2,:); yy=zeros(1,length(Vk1));                                 
Gamma = zeros(1,length(eps)); Lambda = zeros(1,length(eps)); 

% --
% Loop to create the Lambda Functions and Gamma functions. This are
% created by summing each individual Principal Value and ever single
% Delta Function
for i = 1:length(Vk1) 
    % --
    % Creation of the Gaussian

    y = exp(-(eps-Ds(i,i)).^2./(nM^2))./(nM*sqrt(pi));    

    % --
    % Normalization of the Gaussian. Make sure that the integral of the
    % Gaussian is equal to the unity
    y = y./trapz(eps,y); 
    % --
    % yy function is a function that is not incorporated to the final
    % Gamma in case it gives nan. The nan values arises when the trapz
    % function evaluates zero i.e, the step (separation) of eps is not
    % small enough
    yy = 2 * pi * (Vk1(i) * T1A + Vk2(i) * T1B)^2 .* y;

    lambda = 1./(eps - Ds(i,i));
    lambda  = lambda./(max(abs(lambda)));
    Lambda = Lambda + (Vk1(i) * T1A + Vk2(i) * T1B)^2.*lambda;

    % --
    % A check point to make sure we erase any nan data

    if isnan(yy) == 0
        Gamma = Gamma + yy;
    elseif isnan(yy) == 1
        chck(l) = i;
        l = l+1;
    end
end
% --
%Fermi Dirac Distribution:
FDd = 1./(1+exp((eps-ef).*beta));
frct1 = zeros(1,length(ed));
for j=1:1:length(ed)
    % -- 
    % Friction Tensor gamma that will be integrated
    gamma = exp(beta.*(eps-ef)).*(Gamma.*FDd./((eps - ed(j)- Lambda).^2 + (Gamma./2).^2)).^2;
    % Erase all the NaN data from the equation above 
    TF = isnan(gamma);
    gamma(TF) = 0;
    frct1(j) = (beta*g^2)/(2*pi)*trapz(eps,gamma);
end
figure; hold on;
plot(ed,frct1)
ylabel('$\frac{\gamma}{m \omega}$','Interpreter','latex', 'FontSize',12,...
   'FontWeight','bold')
xlabel('Ad-particle energy $\varepsilon_d(x)$','Interpreter','latex')
xlim([edi edf])


%% Friction Tensor in the edges: T parameter
clear variables
tic
beta = 1/0.015; e0 = 0.15; v = 10; w = 14; g = 0.08;
ef = 0; N=800; T1A=linspace(0.085,0.115,200); T1B = T1A./3; psi=zeros(2*N);
% --
%Construction of the unperturbed Hamiltonian:
psi(1,2)=v; psi(2*N,2*N-1)=v;
for ii=2:2:2*(N-1)
    psi(ii,ii-1)=v;
    psi(ii,ii+1)=w;
    psi(ii+1,ii)=w;
    psi(ii+1,ii+2)=v;
end
% --
% Diagonalization of the matrix
[V,D]=eig(psi);
% --
% Sort of the eigenstates and energies
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
% --
% Linewidth of the Deltas function will be between the finite
% difference between the two states with closest energy difference: the
% first and second state
clear y Gamma Vk; nM=0.4*(Ds(2,2)-Ds(1,1)); l = 1; chck = 0;
% --
% Number of steps taken in consideration for the construction of the
% Gamma and Lambda Function. This term is the more sensitive and a small
% change gives unsense results


% --
% Time steps taken to evaluate the integration
ed = 0.01;
eps = linspace(Ds(1,1),Ds(2*N,2*N),(Ds(2*N,2*N)-Ds(1,1))/(4*(Ds(2,2)-Ds(1,1)))); 
% --
% The term V(1,:) gives the projection of |1_A> into the d_k^\dagger. 
% The other term V(:,1) gives the eigenstate with energy E_1

Vk1 = Vs(1,:); Vk2 = Vs(2,:); yy=zeros(1,length(Vk1));                                 

% --
% Loop to create the Lambda Functions and Gamma functions. This are
% created by summing each individual Principal Value and ever single
% Delta Function
frct = zeros(1,length(T1A));
for ii = 1:length(T1A)
    Gamma = zeros(1,length(eps)); Lambda = zeros(1,length(eps)); 
    for i = 1:length(Vk1) 
        % --
        % Creation of the Gaussian
        y = exp(-(eps-Ds(i,i)).^2./(nM^2))./(nM*sqrt(pi));    
    
        % --
        % Normalization of the Gaussian. Make sure that the integral of the
        % Gaussian is equal to the unity
        y = y./trapz(eps,y); 
        % --
        % yy function is a function that is not incorporated to the final
        % Gamma in case it gives nan. The nan values arises when the trapz
        % function evaluates zero i.e, the step (separation) of eps is not
        % small enough
        yy = 2 * pi * (Vk1(i) * T1A(ii) + Vk2(i) * T1B(ii))^2 .* y;
    
        lambda = 1./(eps - Ds(i,i));
        lambda  = lambda./(max(abs(lambda)));
        Lambda = Lambda + (Vk1(i) * T1A(ii) + Vk2(i) * T1B(ii))^2.*lambda;
    
        % --
        % A check point to make sure we erase any nan data
    
        if isnan(yy) == 0
            Gamma = Gamma + yy;
        elseif isnan(yy) == 1
            chck(l) = i;
            l = l+1;
        end
    end
    % --
    %Fermi Dirac Distribution:
    FDd = 1./(1+exp((eps-ef).*beta));
    % -- 
    % Friction Tensor gamma that will be integrated
    gamma = exp(beta.*(eps-ef)).*(Gamma.*FDd./((eps - ed- Lambda).^2 + (Gamma./2).^2)).^2;
    % Erase all the NaN data from the equation above 
    TF = isnan(gamma);
    gamma(TF) = 0;
    frct(ii) = (beta*g^2)/(2*pi)*trapz(eps,gamma);
end
figure; hold on;
plot(T1A,frct)
ylabel('$\frac{\gamma}{m \omega}$','Interpreter','latex', 'FontSize',12,...
   'FontWeight','bold')
xlabel('Ad-particle energy $\varepsilon_d(x)$','Interpreter','latex')

%% Friction Tensor in the Bulk: Metallic phase
clear variables
tic
beta = 1/0.015; e0 = 0.15; v = 10; w = 10; g = 0.08; 
ef = 0; N=800; TmA=0.1; TmB = TmA/3; psi=zeros(2*N); 
edi = e0-0.45; edf = e0+0.15;
% --
%Construction of the unperturbed Hamiltonian:
psi(1,2)=v; psi(2*N,2*N-1)=v;
for ii=2:2:2*(N-1)
    psi(ii,ii-1)=v;
    psi(ii,ii+1)=w;
    psi(ii+1,ii)=w;
    psi(ii+1,ii+2)=v;
end
% --
% Diagonalization of the matrix
[V,D]=eig(psi);
% --
% Sort of the eigenstates and energies
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
% --
% Linewidth of the Deltas function will be between the finite
% difference between the two states with closest energy difference: the
% first and second state
clear y Gamma Vk; nM=0.4*(Ds(2,2)-Ds(1,1)); l = 1; chck = 0;

% --
% Number of steps taken in consideration for the construction of the
% Gamma and Lambda Function. This term is the more sensitive and a small
% change gives unsense results


% --
% Time steps taken to evaluate the integration
ed = g*sqrt(2).*linspace((edi-e0)/(sqrt(2)*g),(edf-e0)/(sqrt(2)*g),(edf-edi)/(1*(Ds(2,2)-Ds(1,1)))) + e0;
eps = linspace(Ds(1,1),Ds(2*N,2*N),(Ds(2*N,2*N)-Ds(1,1))/(0.5*(Ds(2,2)-Ds(1,1)))); 

% --
% The term V(1,:) gives the projection of |1_A> into the d_k^\dagger. 
% The other term V(:,1) gives the eigenstate with energy E_1

VkA = Vs(N+1,:); VkB = Vs(N,:); VkBB = Vs(N+2,:); yy=zeros(1,length(VkA));
Gamma = zeros(1,length(eps)); Lambda = zeros(1,length(eps)); 

% --
% Loop to create the Lambda Functions and Gamma functions. This are
% created by summing each individual Principal Value and ever single
% Delta Function
for i = 1:length(VkA) 
    % --
    % Creation of the Gaussian
    y = exp(-(eps-Ds(i,i)).^2./(nM^2))./(nM*sqrt(pi)); 

    % --
    % Normalization of the Gaussian. Make sure that the integral of the
    % Gaussian is equal to the unity
    y = y./trapz(eps,y); 
    % --
    % yy function is a function that is not incorporated to the final
    % Gamma in case it gives nan. The nan values arises when the trapz
    % function evaluates zero i.e, the step (separation) of eps is not
    % small enough
    yy = 2 * pi * (VkA(i) * TmA + VkB(i) * TmB + VkBB(i) * TmB)^2 .* y;


    lambda = 1./(eps - Ds(i,i));
    lambda  = lambda./(max(abs(lambda)));
    Lambda = Lambda + (VkA(i) * TmA + VkB(i) * TmB + VkBB(i) * TmB)^2.*lambda;


    % --
    % A check point to make sure we erase any nan data

    if isnan(yy) == 0
        Gamma = Gamma + yy;
    elseif isnan(yy) == 1
        chck(l) = i;
        l = l+1;
    end
end
% --
%Fermi Dirac Distribution:
FDd = 1./(1+exp((eps-ef).*beta));
frct = zeros(1,length(ed));
for j=1:1:length(ed)
    % -- 
    % Friction Tensor gamma that will be integrated
    gamma = exp(beta.*(eps-ef)).*(Gamma.*FDd./((eps - ed(j)- Lambda).^2 + (Gamma./2).^2)).^2;
    % Erase all the NaN data from the equation above 
    TF = isnan(gamma);
    gamma(TF) = 0;
    frct(j) = (beta*g^2)/(2*pi)*trapz(eps,gamma);
end
figure; hold on;
plot(ed,frct)
ylabel('$\frac{\gamma}{m \omega}$','Interpreter','latex', 'FontSize',12,...
   'FontWeight','bold')
xlabel('Ad-particle energy $\varepsilon_0(R,x_B)$','Interpreter','latex')
xlim([edi edf])

%% Double axis plot
figure; hold on;
yyaxis left
plot(ed,frct)
ylabel('$\frac{\gamma}{m \omega}$','Interpreter','latex', 'FontSize',12,...
   'FontWeight','bold')
xlabel('Ad-particle energy $\varepsilon_d(x)$','Interpreter','latex')
xlim([edi edf])
plot(ed((length(frct)+1)/2),frct((length(frct)+1)/2),"r*")

yyaxis right
plot(ed,frct1)
ylabel('$\frac{\gamma}{m \omega}$','Interpreter','latex', 'FontSize',12,...
   'FontWeight','bold')
xlabel('Ad-particle energy $\varepsilon_d(x)$','Interpreter','latex')
xlim([edi edf])
plot(ed((length(frct1)+1)/2),frct1((length(frct1)+1)/2),"r*")
