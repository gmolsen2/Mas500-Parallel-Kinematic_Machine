close all;
clear;
clc;

%Test conditions
n_i = 20; %number of discretizations for z-direction
n_j = 20; %number of discretizations for y-direction

px=0;  %[deg]
py=80; %[deg]
pz=0;  %[deg]
A = A_rot(px,py,pz);

%Importing dimensions from external file
Constants = const();
A1m = Constants(1:3);
B1m = Constants(4:6);
C1m = Constants(7:9);
D1m = Constants(10:12);
E1m = Constants(13:15);
F1m = Constants(16:18);
A2_y = Constants(19); 
A2_z = Constants(20);
B2_y = Constants(21);
B2_z = Constants(22);
C2_y = Constants(23);
C2_z = Constants(24);
D2_y = Constants(25);
D2_z = Constants(26);
E2_y = Constants(27);
E2_z = Constants(28);
F2_y = Constants(29);
F2_z = Constants(30);
CoR_Ry = Constants(37);
CoR_Rz = Constants(38);
CoR_Ly = Constants(39); 
CoR_Lz = Constants(40);

%Linkage geometry
L_long = 1.61; %[m]
L_short = 1.45; %[m]
k_long = 78.2; %[N/mu m]
k_short = 81.2; %[N/mu m]

%Additional frame geometry
S_C2m = [0.23252354; -0.043; 0.125]; %[m]
S_D2m = [0.14752354; -0.12; -0.125]; %[m]
S_E2m = [0.14752354; 0.12; -0.125]; %[m]
S_F2m = [0.23252354; 0.043; 0.125]; %[m]

%For loop parameters
AxesY = linspace(-500,500,n_i);
AxesZ = linspace(0,920,n_j);
Resultsx = zeros(n_i,n_j);
Resultsy = zeros(n_i,n_j);
Resultsz = zeros(n_i,n_j);
SingularityCheck = zeros(n_i,n_j);


%Checking stiffness of new system against old values
for i=1:n_i
    z = AxesZ(i) * 1E-3;
    for j=1:n_j
        y = AxesY(j) * 1E-3;
        TCP = [0; y; z];

        %Calculations
        %Finding global position of ball socket joint on tool platform
        A1 = TCP + A*A1m;
        B1 = TCP + A*B1m;
        C1 = TCP + A*C1m;
        D1 = TCP + A*D1m;
        E1 = TCP + A*E1m;
        F1 = TCP + A*F1m;
        
        %Finding vectors from CoR to joints in global coordinate system
        S_C2 = A*S_C2m;
        S_D2 = A*S_D2m;
        S_E2 = A*S_E2m;
        S_F2 = A*S_F2m;
        
        %Finding x-position of rotational point on carriage
        A2_x = -sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2) + A1(1);
        B2_x = A2_x;
        CoR_Rx = -sqrt(L_short^2 - (CoR_Ry+S_D2(2)-D1(2))^2 - (CoR_Rz+S_D2(3)-D1(3))^2) - S_D2(1) + D1(1);
        CoR_Lx = -sqrt(L_short^2 - (CoR_Ly+S_E2(2)-E1(2))^2 - (CoR_Lz+S_E2(3)-E1(3))^2) - S_E2(1) + E1(1);
        
        %Resulting ficture positions
        CoR_R = [CoR_Rx; CoR_Ry; CoR_Rz];
        CoR_L = [CoR_Lx; CoR_Ly; CoR_Lz];
        A2 = [A2_x; A2_y; A2_z];
        B2 = [B2_x; B2_y; B2_z];
        C2 = CoR_R + S_C2;
        D2 = CoR_R + S_D2;
        E2 = CoR_L + S_E2;
        F2 = CoR_L + S_F2;
        
        %Finding directional vectors of arms
        Ua = (A2-A1)/norm(A2-A1);
        Ub = (B2-B1)/norm(B2-B1);
        Uc = (C2-C1)/norm(C2-C1);
        Ud = (D2-D1)/norm(D2-D1);
        Ue = (E2-E1)/norm(E2-E1);
        Uf = (F2-F1)/norm(F2-F1);

        %Finding moment arm at ball socket joints from TCP
        rA = cross((A1-TCP),Ua);
        rB = cross((B1-TCP),Ub);
        rC = cross((C1-TCP),Uc);
        rD = cross((D1-TCP),Ud);
        rE = cross((E1-TCP),Ue);
        rF = cross((F1-TCP),Uf);

        %Static matrix H
        H = [Ua, Ub, Uc, Ud Ue, Uf; rA, rB, rC, rD, rE, rF];

        %Cartesian stiffness matrix
        k_L = [k_long, k_long, k_short, k_short, k_short, k_short];
        k_Cart = H*diag(k_L)*H';

        %Deflection of rods due to moment pr. force at TCP
        delta_x = k_Cart \ [1;0;0;0;0;0];
        delta_y = k_Cart \ [0;1;0;0;0;0];
        delta_z = k_Cart \ [0;0;1;0;0;0];

        %Stiffness in global cartesian direction
        k_x = 1/norm(delta_x(1:3));
        k_y = 1/norm(delta_y(1:3));
        k_z = 1/norm(delta_z(1:3));
        
        if isreal(CoR_Rx+CoR_Lx)
            Resultsx(i,j) = k_x;
            Resultsy(i,j) = k_y;
            Resultsz(i,j) = k_z;
            SingularityCheck(i,j) = cond(H);
        end
    end
end

figure(1)
surf(AxesY,AxesZ,Resultsx)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Stiffness [N/\mu m]')
title('Stiffness in X-direction')

figure(2)
surf(AxesY,AxesZ,Resultsy)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Stiffness [N/\mu m]')
title('Stiffness in Y-direction')

figure(3)
surf(AxesY,AxesZ,Resultsz)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Stiffness [N/\mu m]')
title('Stiffness in Z-direction')

figure(4)
surf(AxesY,AxesZ,SingularityCheck)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Condition number')
title('Condition number for matrix inversion')

figure(200)
plot3(A1(1),A1(2),A1(3),'.' ,A2(1), A2(2), A2(3), '.','MarkerSize',15, 'color','k')
hold on
xlabel('X-axis [m]')
ylabel('Y-axis [m]')
zlabel('Z-axis [m]')
text(A1(1)-0.05,A1(2)+0.1,A1(3)+0.05,'A1', 'FontSize',8)
text(A2(1)+0.05,A2(2)+0.05,A2(3)+0.05,'A2', 'FontSize',8)
plot3(B1(1),B1(2),B1(3),'.' ,B2(1), B2(2), B2(3), '.','MarkerSize',15, 'color','k')
text(B1(1)-0.05,B1(2)+0.1,B1(3)+0.05,'B1', 'FontSize',8)
text(B2(1)+0.05,B2(2)+0.05,B2(3)+0.05,'B2', 'FontSize',8)
plot3(C1(1),C1(2),C1(3),'.' ,C2(1), C2(2), C2(3), '.','MarkerSize',15, 'color','k')
text(C1(1)-0.05,C1(2)+0.1,C1(3)+0.05,'C1', 'FontSize',8)
text(C2(1)+0.05,C2(2)+0.05,C2(3)+0.05,'C2', 'FontSize',8)
plot3(D1(1),D1(2),D1(3),'.' ,D2(1), D2(2), D2(3), '.','MarkerSize',15, 'color','k')
text(D1(1)-0.05,D1(2)+0.1,D1(3)+0.05,'D1', 'FontSize',8)
text(D2(1)+0.05,D2(2)+0.05,D2(3)+0.05,'D2', 'FontSize',8)
plot3(E1(1),E1(2),E1(3),'.' ,E2(1), E2(2), E2(3), '.','MarkerSize',15, 'color','k')
text(E1(1)-0.05,E1(2)+0.1,E1(3)+0.05,'E1', 'FontSize',8)
text(E2(1)+0.05,E2(2)+0.05,E2(3)+0.05,'E2', 'FontSize',8)
plot3(F1(1),F1(2),F1(3),'.' ,F2(1), F2(2), F2(3), '.','MarkerSize',15, 'color','k')
text(F1(1)-0.05,F1(2)+0.1,F1(3)+0.05,'F1', 'FontSize',8)
text(F2(1)+0.05,F2(2)+0.05,F2(3)+0.05,'F2', 'FontSize',8)
plot3(TCP(1),TCP(2),TCP(3), '.','MarkerSize',15, 'color','k')
text(TCP(1)-0.05,TCP(2)+0.15,TCP(3)-0.05,'TCP', 'FontSize',8)
plot3(CoR_R(1),CoR_R(2),CoR_R(3),'.' ,CoR_L(1), CoR_L(2), CoR_L(3), '.','MarkerSize',15, 'color','k')
text(CoR_R(1)-0.05,CoR_R(2)+0.1,CoR_R(3)+0.05,'CoR_R', 'FontSize',8)
text(CoR_L(1)+0.05,CoR_L(2)+0.05,CoR_L(3)+0.05,'CoR_L', 'FontSize',8)

line([A1(1), A2(1)], [A1(2), A2(2)],[A1(3), A2(3)], 'Linewidth', 2, 'color', '#0072BD')
line([B1(1), B2(1)], [B1(2), B2(2)],[B1(3), B2(3)], 'Linewidth', 2, 'color', '#D95319')
line([C1(1), C2(1)], [C1(2), C2(2)],[C1(3), C2(3)], 'Linewidth', 2, 'color', '#EDB120')
line([D1(1), D2(1)], [D1(2), D2(2)],[D1(3), D2(3)], 'Linewidth', 2, 'color', '#7E2F8E')
line([E1(1), E2(1)], [E1(2), E2(2)],[E1(3), E2(3)], 'Linewidth', 2, 'color', '#77AC30')
line([F1(1), F2(1)], [F1(2), F2(2)],[F1(3), F2(3)], 'Linewidth', 2, 'color', '#4DBEEE')
TCP2 = TCP + A*[0;0;0.3];
line([TCP(1), TCP2(1)], [TCP(2), TCP2(2)],[TCP(3), TCP2(3)], 'Linewidth', 2, 'color', 'k')

function y = A_rot(px,py,pz)
    %Rotation matrix. Bryant angles chosen is [x-y-z]
    Ax = [1, 0, 0;
          0, cosd(px), -sind(px);
          0, sind(px), cosd(px)];
    Ay = [cosd(py), 0, sind(py);
          0, 1, 0;
          -sind(py), 0, cosd(py)];
    Az = [cosd(pz), -sind(pz), 0;
          sind(pz), cosd(pz), 0;
          0, 0, 1];
    
    y = Ax * Ay * Az;

end
