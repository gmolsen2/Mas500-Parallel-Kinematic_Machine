close all;
clear;
clc;

%Test conditions
n_i = 20; %number of discretizations for z-direction
n_j = 20; %number of discretizations for y-direction

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
P1m = Constants(31:33);
P2m = Constants(34:36);
d_P = norm(P2m-P1m);

%Linkage geometry
L_long = 1.61; %[m]
L_short = 1.45; %[m]
k_long = 78.2; %[N/mu m]
k_short = 81.2; %[N/mu m]
d = 0.079; %[m]

%Matrix creation
Resultsx_pos = zeros(n_i,n_j);
Resultsx_neg = zeros(n_i,n_j);
Extensionx_posC = zeros(n_i,n_j);
Extensionx_posF = zeros(n_i,n_j);
Extensionx_negC = zeros(n_i,n_j);
Extensionx_negF = zeros(n_i,n_j);
AxesY = linspace(-500,500,n_i);
AxesZ = linspace(0,920,n_j);

%Checking stiffness of new system against old values
for i=1:n_i
    z = AxesZ(i) * 1E-3;
    for j=1:n_j
        y = AxesY(j) * 1E-3;
        TCP = [0; y; z];

        %Calculating positive X
        px = 0; breaksignal = boolean(0);

        while px < 90 && breaksignal == boolean(0)
            px = px + 1; error = 1; count = 0;

            %Initial guess
            pz = -px;

            %NR-solver finding z rotation of head and x position of top gantry
            while error > 1E-5 && count < 1000 && breaksignal == boolean(0)
                A = A_rot(px,0,pz);
                A1 = TCP + A*A1m;
                B1 = TCP + A*B1m;

                dA1xdpz = -A1m(1)*sind(pz) - A1m(2)*cosd(pz);
                dA1ydpz = A1m(1)*cosd(px)*cosd(pz) - A1m(2)*cosd(px)*sind(pz);
                dA1zdpz = A1m(1)*sind(px)*cosd(pz) - A1m(2)*sind(px)*sind(pz);

                dB1xdpz = -B1m(1)*sind(pz) - B1m(2)*cosd(pz);
                dB1ydpz = B1m(1)*cosd(px)*cosd(pz) - B1m(2)*cosd(px)*sind(pz);
                dB1zdpz = B1m(1)*sind(px)*cosd(pz) - B1m(2)*sind(px)*sind(pz);

                f1 = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2 = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f = f1 - f2;

                df1 = dA1ydpz * (A2_y-A1(2)) + dA1zdpz * (A2_z-A1(3));
                df2 = sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                df3 = dB1ydpz * (B2_y-B1(2)) + dB1zdpz * (B2_z-B1(3));
                df4 = sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                df = dA1xdpz - df1/df2 - dB1xdpz + df3/df4;

                pz = pz - f / df;

                A = A_rot(px,0,pz);
                A1 = TCP + A*A1m;
                B1 = TCP + A*B1m;
                f1_new = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2_new = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f_new = f1_new - f2_new;

                count = count + 1;
                error = norm(f_new);

                if count > 50 && error > 0.2
                    breaksignal = boolean(1);
                end
            end

            A = A_rot(px,0,pz);

            A1 = TCP + A*A1m;
            B1 = TCP + A*B1m;
            C1 = TCP + A*C1m;
            D1 = TCP + A*D1m;
            E1 = TCP + A*E1m;
            F1 = TCP + A*F1m;
            
            %Finding x-position of ball socet joint on carriage
            A2_x = -sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2) + A1(1);
            B2_x = A2_x;
            D2_x = -sqrt(L_short^2 - (D2_y-D1(2))^2 - (D2_z-D1(3))^2) + D1(1);
            C2_x = D2_x + 0.0085;
            E2_x = -sqrt(L_short^2 - (E2_y-E1(2))^2 - (E2_z-E1(3))^2) + E1(1);
            F2_x = E2_x + 0.0085;
            Extention_C = sqrt((C2_x-C1(1))^2 + (C2_y-C1(2))^2 + (C2_z-C1(3))^2) - L_short;
            Extention_F = sqrt((F2_x-F1(1))^2 + (F2_y-F1(2))^2 + (F2_z-F1(3))^2) - L_short;
            
            %Resulting fixture position
            A2 = [A2_x; A2_y; A2_z];
            B2 = [B2_x; B2_y; B2_z];
            C2 = [C2_x; C2_y; C2_z];
            D2 = [D2_x; D2_y; D2_z];
            E2 = [E2_x; E2_y; E2_z];
            F2 = [F2_x; F2_y; F2_z];

            %Finding distances between vectors
            v_A = A2 - A1;
            v_B = B2 - B1;
            v_C = C2 - C1;
            v_F = F2 - F1;

            d_AF = D_lineline(A1,v_A,F1,v_F);
            d_BC = D_lineline(B1,v_B,C1,v_C);

            %Checking for irrelevant results
            if ~isreal(pz)
                break
            elseif ~isreal(D2_x+E2_x+Extention_C+Extention_F)
                break
            elseif d_AF < d | d_BC < d
                break
            end
            
            %Saving results if NR solver is within error
            if breaksignal == boolean(0)
                Resultsx_pos(i,j) = px;
                Extensionx_posC(i,j) = Extention_C * 1E3;
                Extensionx_posF(i,j) = Extention_F * 1E3;
            end

        end

        %Calculating negative X
        px = 0; breaksignal = boolean(0);

        while px > -90 && breaksignal == boolean(0)
            px = px - 1; error = 1; count = 0;

            %Initial guess
            pz = -px;

            %NR-solver finding z rotation of head and x position of top gantry
            while error > 1E-5 && count < 1000 && breaksignal == boolean(0)
                A = A_rot(px,0,pz);
                A1 = TCP + A*A1m;
                B1 = TCP + A*B1m;

                dA1xdpz = -A1m(1)*sind(pz) - A1m(2)*cosd(pz);
                dA1ydpz = A1m(1)*cosd(px)*cosd(pz) - A1m(2)*cosd(px)*sind(pz);
                dA1zdpz = A1m(1)*sind(px)*cosd(pz) - A1m(2)*sind(px)*sind(pz);

                dB1xdpz = -B1m(1)*sind(pz) - B1m(2)*cosd(pz);
                dB1ydpz = B1m(1)*cosd(px)*cosd(pz) - B1m(2)*cosd(px)*sind(pz);
                dB1zdpz = B1m(1)*sind(px)*cosd(pz) - B1m(2)*sind(px)*sind(pz);

                f1 = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2 = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f = f1 - f2;

                df1 = dA1ydpz * (A2_y-A1(2)) + dA1zdpz * (A2_z-A1(3));
                df2 = sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                df3 = dB1ydpz * (B2_y-B1(2)) + dB1zdpz * (B2_z-B1(3));
                df4 = sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                df = dA1xdpz - df1/df2 - dB1xdpz + df3/df4;

                pz = pz - f / df;

                A = A_rot(px,0,pz);
                A1 = TCP + A*A1m;
                B1 = TCP + A*B1m;
                f1_new = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2_new = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f_new = f1_new - f2_new;

                count = count + 1;
                error = norm(f_new);

                if count > 50 && error > 0.2
                    breaksignal = boolean(1);
                end
            end

            A = A_rot(px,0,pz);

            A1 = TCP + A*A1m;
            B1 = TCP + A*B1m;
            C1 = TCP + A*C1m;
            D1 = TCP + A*D1m;
            E1 = TCP + A*E1m;
            F1 = TCP + A*F1m;
            
            %Finding x-position of ball socet joint on carriage
            A2_x = -sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2) + A1(1);
            B2_x = A2_x;
            D2_x = -sqrt(L_short^2 - (D2_y-D1(2))^2 - (D2_z-D1(3))^2) + D1(1);
            C2_x = D2_x + 0.0085;
            E2_x = -sqrt(L_short^2 - (E2_y-E1(2))^2 - (E2_z-E1(3))^2) + E1(1);
            F2_x = E2_x + 0.0085;
            Extention_C = sqrt((C2_x-C1(1))^2 + (C2_y-C1(2))^2 + (C2_z-C1(3))^2) - L_short;
            Extention_F = sqrt((F2_x-F1(1))^2 + (F2_y-F1(2))^2 + (F2_z-F1(3))^2) - L_short;
            
            %Resulting fixture position
            A2 = [A2_x; A2_y; A2_z];
            B2 = [B2_x; B2_y; B2_z];
            C2 = [C2_x; C2_y; C2_z];
            D2 = [D2_x; D2_y; D2_z];
            E2 = [E2_x; E2_y; E2_z];
            F2 = [F2_x; F2_y; F2_z];

            %Finding distances between vectors
            v_A = A2 - A1;
            v_B = B2 - B1;
            v_C = C2 - C1;
            v_F = F2 - F1;

            d_AF = D_lineline(A1,v_A,F1,v_F);
            d_BC = D_lineline(B1,v_B,C1,v_C);

            %Checking for irrelevant results
            if ~isreal(pz)
                break
            elseif ~isreal(D2_x+E2_x+Extention_C+Extention_F)
                break
            elseif d_AF < d | d_BC < d
                break
            end

            %Saving results if NR solver is within error
            if breaksignal == boolean(0)
                Resultsx_neg(i,j) = -px;
                Extensionx_negC(i,j) = Extention_C * 1E3;
                Extensionx_negF(i,j) = Extention_F * 1E3;
            end

        end
    end
end

figure(1)
surf(AxesY,AxesZ,Resultsx_pos)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('rotation [deg]')
title('Max rotation positive X')

% figure(2)
% surf(AxesY,AxesZ,Extensionx_posC)
% grid on
% colorbar
% xlabel('Y-axis [mm]')
% ylabel('Z-axis [mm]')
% zlabel('Joint coordinate C [mm]')
% title('Maximum Extention C')
% 
% figure(3)
% surf(AxesY,AxesZ,Extensionx_posF)
% grid on
% colorbar
% xlabel('Y-axis [mm]')
% ylabel('Z-axis [mm]')
% zlabel('Joint coordinate F [mm]')
% title('Maximum Extention F')

figure(4)
surf(AxesY,AxesZ,Resultsx_neg)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('rotation [deg]')
title('Max rotation negative X')

% figure(5)
% surf(AxesY,AxesZ,Extensionx_negC)
% grid on
% colorbar
% xlabel('Y-axis [mm]')
% ylabel('Z-axis [mm]')
% zlabel('Joint coordinate C [mm]')
% title('Maximum Extention C')
% 
% figure(6)
% surf(AxesY,AxesZ,Extensionx_negF)
% grid on
% colorbar
% xlabel('Y-axis [mm]')
% ylabel('Z-axis [mm]')
% zlabel('Joint coordinate F [mm]')
% title('Maximum Extention F')

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

function y = D_lineline(point1,vector1,point2,vector2)
    %This function computes a distance value between two non paralell lines
    cr = cross(vector1,vector2);
    pointpoint = point2 - point1;
    y = abs(dot(cr/norm(cr),pointpoint));
end

function y = D_linepoint(point,pointline,vectorline)
    %This function computes a distance value between a point and a line
    pointpoint = pointline - point;
    y = norm(cross(pointpoint,vectorline) / norm(vectorline));
end