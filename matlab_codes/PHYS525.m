clear all; 
t=1;
kx=linspace(-2*pi,2*pi,20); ky=linspace(-2*pi,2*pi,20);
for jj=1:length(kx)
    for j=1:length(ky)
        m=-t*[2*cos(ky(j)) 1. exp(-1i*kx(jj));1. 2*cos(ky(j)+2*pi/3) 1;exp(1i*kx(jj)) 1 2*cos(ky(j)+4*pi/3)];
        [V,D]=eig(m);
        band1(j,jj)=D(1,1);
        band2(j,jj)=D(2,2);
        band3(j,jj)=D(3,3);
    end
end
[KX,KY]=meshgrid(kx,ky);

%%
clear all; 
t=1; qq=30; rr=1; tt=1; list=0;
kx=linspace(-2*pi,2*pi,30); ky=linspace(-2*pi,2*pi,30); figure; hold on
for q=3:qq
    for p=q+1:qq+1
        if ismember(q/p,list) == 0
            list(tt)=q/p; tt=tt+1;
            for jj=1:length(kx)
                psi=zeros(q);
                psi(1,q)=-t*exp(1i*kx(jj)); psi(q,1)=-t*exp(-1i*kx(jj));
                psi(1,2)=-t*(psi(1,2)+1);  psi(q,q-1)=-t*(psi(q,q-1)+1);
                for ii=1:length(ky)
                    psi(1,1)=-2*t*cos(ky(ii)); 
                    psi(q,q)=-2*t*cos(ky(ii)+2*pi*(q-1)*p/q);
                    for ll=2:1:q-1
                        psi(ll,ll)=-2*t*cos(ky(ii)+2*(ll-1)*pi/q);
                        psi(ll,ll-1)=-t; psi(ll,ll+1)=-t;
                    end
                    [V,D]=eig(psi);
                    for ll=1:length(D)
                        fnloop(1,rr)=D(ll,ll);
                        fnloop(2,rr)=q/p;
                        rr=rr+1;
                    end
                end
            end
        end
    end
end
scatter(fnloop(2,:),fnloop(1,:),0.5,"filled",'k')
hold on
scatter(1-fnloop(2,:),fnloop(1,:),0.5,"filled",'k')


%%
x=linspace(1,10,10); y=linspace(2,10,5);
for i=1:length(x)
    if ismember(x(i),y)==1
        x(i)
    end
end