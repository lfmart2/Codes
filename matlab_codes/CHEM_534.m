clear all; clc;
theta = deg2rad(104.5); r = 0.958;
m_D = 2.0141; m_H = 1.00797; m_O = 15.994;
M = [m_D m_H m_O];
r_D = [-r*sin(theta/2)*(1-m_H/(m_O+3*m_H)) -r*cos(theta/2)*(1-3*m_H/(m_O+3*m_H)) 0];
r_H = [r*sin(theta/2)*(1+m_H/(m_O+3*m_H)) -r*cos(theta/2)*(1-3*m_H/(m_O+3*m_H)) 0];
r_O = [m_H*r*sin(theta/2)/(m_O + 3*m_H) 3*m_H*r*cos(theta/2)/(m_O + 3*m_H) 0];
R = [r_D;r_H;r_O];
I_x = 0; I_y = 0; I_z = 0;
I_xy = 0; I_xz = 0; I_yz = 0;
for i=1:length(M)
    I_x = I_x + M(i)*(R(i,2)^2+R(i,3)^2);
    I_y = I_y + M(i)*(R(i,1)^2+R(i,3)^2);
    I_z = I_z + M(i)*(R(i,1)^2+R(i,2)^2);
    I_xy = I_xy + M(i)*R(i,1)*R(i,2);
    I_xz = I_xz + M(i)*R(i,1)*R(i,3);
    I_yz = I_yz + M(i)*R(i,2)*R(i,3);
end
I = [I_x -I_xy -I_xz;-I_xy I_y I_yz;-I_xz -I_yz I_z];
[V,D]=eig(I);

atm = 1.6605402E-27;
h = 6.62607015e-34;
angs = 1e-10;
A= h/(D(1,1)*atm*angs^2*8*pi^2)
B= h/(D(2,2)*atm*angs^2*8*pi^2)
C= h/(D(3,3)*atm*angs^2*8*pi^2)