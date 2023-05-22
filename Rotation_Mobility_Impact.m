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
CoR_Ry = Constants(37);
CoR_Rz = Constants(38);
CoR_Ly = Constants(39); 
CoR_Lz = Constants(40);
d_P = norm(P2m-P1m);

%Linkage geometry
L_long = 1.61; %[m]
L_short = 1.45; %[m]
k_long = 78.2; %[N/mu m]
k_short = 81.2; %[N/mu m]
d = 0.079; %[m]
r = 0.2; %[m]

%Extra parameters
S_C2m = [0.23252354; -0.043; 0.125]; %[m]
S_D2m = [0.14752354; -0.12; -0.125]; %[m]
S_E2m = [0.14752354; 0.12; -0.125]; %[m]
S_F2m = [0.23252354; 0.043; 0.125]; %[m]

%Matrix creation
Resultsx_pos = zeros(n_i,n_j);
Resultsx_neg = zeros(n_i,n_j);
Resultsy_pos = zeros(n_i,n_j);
Resultsy_neg = zeros(n_i,n_j);
JointCoordx_posE = zeros(n_i,n_j);
JointCoordx_posW = zeros(n_i,n_j);
JointCoordx_negE = zeros(n_i,n_j);
JointCoordx_negW = zeros(n_i,n_j);
JointCoordy_posE = zeros(n_i,n_j);
JointCoordy_posW = zeros(n_i,n_j);
JointCoordy_negE = zeros(n_i,n_j);
JointCoordy_negW = zeros(n_i,n_j);
AxesY = linspace(-500,500,n_i);
AxesZ = linspace(0,920,n_j);

%Checking stiffness of new system against old values
for i=1:n_i
    z = AxesZ(i) * 1E-3;
    for j=1:n_j
        y = AxesY(j) * 1E-3;
        TCP = [0; y; z];
 
        %Calculating positive X
        px = 0; py = 0; breaksignal = boolean(0); Test_C = boolean(1);
        Test_D = boolean(1); Test_E = boolean(1); Test_F = boolean(1);

        while px < 90 && breaksignal == boolean(0) && Test_C && Test_D && Test_E && Test_F
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

                f1_new = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2_new = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f_new = f1_new - f2_new;

                count = count + 1;
                error = norm(f_new);

                if count > 50 && error > 0.2
                    breaksignal = boolean(1);
                end
            end

            %Calculations
            %Finding global position of ball socket joint on tool platform
            A = A_rot(px,0,pz);
            A1 = TCP + A*A1m;
            B1 = TCP + A*B1m;
            C1 = TCP + A*C1m;
            D1 = TCP + A*D1m;
            E1 = TCP + A*E1m;
            F1 = TCP + A*F1m;


            %NR solver for east sled
            errorE = 1; CountE = 0; RE = 0;
            while errorE > 1E-5 && CountE < 1E4
                Ay = A_rot(0,RE,0);
                S_C2 = Ay*S_C2m;
                S_D2 = Ay*S_D2m;
            
                f1 = C1(1) - S_C2(1);
                f2 = sqrt(L_short^2 - (CoR_Ry+S_C2(2)-C1(2))^2 - (CoR_Rz+S_C2(3)-C1(3))^2);
                f3 = D1(1) - S_D2(1);
                f4 = sqrt(L_short^2 - (CoR_Ry+S_D2(2)-D1(2))^2 - (CoR_Rz+S_D2(3)-D1(3))^2);
                f = f1 - f2 - f3 + f4;
            
                dSC2xdR = -sind(RE) * S_C2m(1) + cosd(RE) * S_C2m(3);
                dSC2zdR = -cosd(RE) * S_C2m(1) - sind(RE) * S_C2m(3);
                dSD2xdR = -sind(RE) * S_D2m(1) + cosd(RE) * S_D2m(3);
                dSD2zdR = -cosd(RE) * S_D2m(1) - sind(RE) * S_D2m(3);
            
                df1 = dSC2zdR * (CoR_Rz + S_C2(2)-C1(2));
                df2 = dSD2zdR * (CoR_Rz + S_D2(2)-D1(2));
                df = -dSC2xdR - 0.5 * df1 / f2 + dSD2xdR + 0.5 * df2 / f4;
            
                RE = RE - f/df;

                %Ensuring joint angle is within one rotation
                if RE < 0 | RE > 360
                    RE = wrapTo360(real(RE));
                else
                    RE = real(RE);
                end
            
                %Ensuring only one possible solution
                if RE < 180 && RE > 90
                    RE = RE + 180;
                elseif RE < 270 && RE > 180
                    RE = RE - 180;
                end
            
                Ay = A_rot(0,RE,0);
                S_C2 = Ay*S_C2m;
                S_D2 = Ay*S_D2m;
            
                f1_new = C1(1) - S_C2(1);
                f2_new = sqrt(L_short^2 - (CoR_Ry+S_C2(2)-C1(2))^2 - (CoR_Rz+S_C2(3)-C1(3))^2);
                f3_new = D1(1) - S_D2(1);
                f4_new = sqrt(L_short^2 - (CoR_Ry+S_D2(2)-D1(2))^2 - (CoR_Rz+S_D2(3)-D1(3))^2);
                f_new = f1_new - f2_new - f3_new + f4_new;
            
                errorE = abs(f_new);
                CountE = CountE + 1;
            end
            
            
            %NR solver for west sled
            errorW = 1; CountW = 0; RW = 0;
            while errorW > 1E-5 && CountW < 1E4
                Ay = A_rot(0,RW,0);
                S_E2 = Ay*S_E2m;
                S_F2 = Ay*S_F2m;
            
                f1 = E1(1) - S_E2(1);
                f2 = sqrt(L_short^2 - (CoR_Ly+S_E2(2)-E1(2))^2 - (CoR_Lz+S_E2(3)-E1(3))^2);
                f3 = F1(1) - S_F2(1);
                f4 = sqrt(L_short^2 - (CoR_Ly+S_F2(2)-F1(2))^2 - (CoR_Lz+S_F2(3)-F1(3))^2);
                f = f1 - f2 - f3 + f4;
            
                dSE2xdR = -sind(RW) * S_E2m(1) + cosd(RW) * S_E2m(3);
                dSE2zdR = -cosd(RW) * S_E2m(1) - sind(RW) * S_E2m(3);
                dSF2xdR = -sind(RW) * S_F2m(1) + cosd(RW) * S_F2m(3);
                dSF2zdR = -cosd(RW) * S_F2m(1) - sind(RW) * S_F2m(3);
            
                df1 = dSE2zdR * (CoR_Lz + S_E2(2) - E1(2));
                df2 = dSF2zdR * (CoR_Lz + S_F2(2) - F1(2));
                df = -dSE2xdR - 0.5 * df1 / f2 + dSF2xdR + 0.5 * df2 / f4;
            
                RW = RW - f/df;

                %Ensuring joint angle is within one rotation
                if RW < 0 | RW > 360
                    RW = wrapTo360(real(RW));
                else
                    RW = real(RW);
                end
            
                %Ensuring only one possible solution
                if RW < 180 && RW > 90
                    RW = RW + 180;
                elseif RW < 270 && RW > 180
                    RW = RW - 180;
                end
            
                Ay = A_rot(0,RW,0);
                S_E2 = Ay*S_E2m;
                S_F2 = Ay*S_F2m;
            
                f1_new = E1(1) - S_E2(1);
                f2_new = sqrt(L_short^2 - (CoR_Ly+S_E2(2)-E1(2))^2 - (CoR_Lz+S_E2(3)-E1(3))^2);
                f3_new = F1(1) - S_F2(1);
                f4_new = sqrt(L_short^2 - (CoR_Ly+S_F2(2)-F1(2))^2 - (CoR_Lz+S_F2(3)-F1(3))^2);
                f_new = f1_new - f2_new - f3_new + f4_new;
            
                errorW = abs(f_new);
                CountW = CountW + 1;
            end
            %Finding vectors from CoR to joints in global coordinate system
            S_C2 = A_rot(0,RE,0)*S_C2m;
            S_D2 = A_rot(0,RE,0)*S_D2m;
            S_E2 = A_rot(0,RW,0)*S_E2m;
            S_F2 = A_rot(0,RW,0)*S_F2m;
            
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
            
            %Test for length of linkages
            L_C = norm(C2-C1);
            L_D = norm(D2-D1);
            L_E = norm(E2-E1);
            L_F = norm(F2-F1);
            Test_C = SolutionCheck(L_C, L_short, 0.001);
            Test_D = SolutionCheck(L_D, L_short, 0.001);
            Test_E = SolutionCheck(L_E, L_short, 0.001);
            Test_F = SolutionCheck(L_F, L_short, 0.001);

            %Finding distances between vectors
            v_A = A2 - A1;
            v_B = B2 - B1;
            v_C = C2 - C1;
            v_D = D2 - D1;
            v_E = E2 - E1;
            v_F = F2 - F1;

            d_AF = D_lineline(A1,v_A,F1,v_F);
            d_BC = D_lineline(B1,v_B,C1,v_C);
            d_CD = D_lineline(C1,v_C,D1,v_D);
            d_EF = D_lineline(E1,v_E,F1,v_F);

            %Checking for irrelevant results
            if ~isreal(pz)
                break
            elseif d_AF < d | d_BC < d | d_CD < d | d_EF < d
                break
            end

            %Plotting if NR solver is within error
            if breaksignal == boolean(0) && Test_C && Test_D && Test_D && Test_E && Test_F
                Resultsx_pos(i,j) = px;
                JointCoordx_posE(i,j) = RE;
                JointCoordx_posW(i,j) = RW;
            end

        end

        %Calculating negative X
        px = 0; py = 0; breaksignal = boolean(0); Test_C = boolean(1);
        Test_D = boolean(1); Test_E = boolean(1); Test_F = boolean(1);

        while px > -90 && breaksignal == boolean(0) && Test_C && Test_D && Test_D && Test_E && Test_F
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

                f1_new = A1(1) - sqrt(L_long^2 - (A2_y-A1(2))^2 - (A2_z-A1(3))^2);
                f2_new = B1(1) - sqrt(L_long^2 - (B2_y-B1(2))^2 - (B2_z-B1(3))^2);
                f_new = f1_new - f2_new;

                count = count + 1;
                error = norm(f_new);

                if count > 50 && error > 0.2
                    breaksignal = boolean(1);
                end
            end

            %Calculations
            %Finding global position of ball socket joint on tool platform
            A = A_rot(px,0,pz);
            A1 = TCP + A*A1m;
            B1 = TCP + A*B1m;
            C1 = TCP + A*C1m;
            D1 = TCP + A*D1m;
            E1 = TCP + A*E1m;
            F1 = TCP + A*F1m;


            %NR solver for east sled
            errorE = 1; CountE = 0; RE = 0;
            while errorE > 1E-5 && CountE < 1E4
                Ay = A_rot(0,RE,0);
                S_C2 = Ay*S_C2m;
                S_D2 = Ay*S_D2m;
            
                f1 = C1(1) - S_C2(1);
                f2 = sqrt(L_short^2 - (CoR_Ry+S_C2(2)-C1(2))^2 - (CoR_Rz+S_C2(3)-C1(3))^2);
                f3 = D1(1) - S_D2(1);
                f4 = sqrt(L_short^2 - (CoR_Ry+S_D2(2)-D1(2))^2 - (CoR_Rz+S_D2(3)-D1(3))^2);
                f = f1 - f2 - f3 + f4;
            
                dSC2xdR = -sind(RE) * S_C2m(1) + cosd(RE) * S_C2m(3);
                dSC2zdR = -cosd(RE) * S_C2m(1) - sind(RE) * S_C2m(3);
                dSD2xdR = -sind(RE) * S_D2m(1) + cosd(RE) * S_D2m(3);
                dSD2zdR = -cosd(RE) * S_D2m(1) - sind(RE) * S_D2m(3);
            
                df1 = dSC2zdR * (CoR_Rz + S_C2(2)-C1(2));
                df2 = dSD2zdR * (CoR_Rz + S_D2(2)-D1(2));
                df = -dSC2xdR - 0.5 * df1 / f2 + dSD2xdR + 0.5 * df2 / f4;
            
                RE = RE - f/df;

                %Ensuring joint angle is within one rotation
                if RE < 0 | RE > 360
                    RE = wrapTo360(real(RE));
                else
                    RE = real(RE);
                end
            
                %Ensuring only one possible solution
                if RE < 180 && RE > 90
                    RE = RE + 180;
                elseif RE < 270 && RE > 180
                    RE = RE - 180;
                end
            
                Ay = A_rot(0,RE,0);
                S_C2 = Ay*S_C2m;
                S_D2 = Ay*S_D2m;
            
                f1_new = C1(1) - S_C2(1);
                f2_new = sqrt(L_short^2 - (CoR_Ry+S_C2(2)-C1(2))^2 - (CoR_Rz+S_C2(3)-C1(3))^2);
                f3_new = D1(1) - S_D2(1);
                f4_new = sqrt(L_short^2 - (CoR_Ry+S_D2(2)-D1(2))^2 - (CoR_Rz+S_D2(3)-D1(3))^2);
                f_new = f1_new - f2_new - f3_new + f4_new;
            
                errorE = abs(f_new);
                CountE = CountE + 1;
            end
            
            
            %NR solver for west sled
            errorW = 1; CountW = 0; RW = 0;
            while errorW > 1E-5 && CountW < 1E4
                Ay = A_rot(0,RW,0);
                S_E2 = Ay*S_E2m;
                S_F2 = Ay*S_F2m;
            
                f1 = E1(1) - S_E2(1);
                f2 = sqrt(L_short^2 - (CoR_Ly+S_E2(2)-E1(2))^2 - (CoR_Lz+S_E2(3)-E1(3))^2);
                f3 = F1(1) - S_F2(1);
                f4 = sqrt(L_short^2 - (CoR_Ly+S_F2(2)-F1(2))^2 - (CoR_Lz+S_F2(3)-F1(3))^2);
                f = f1 - f2 - f3 + f4;
            
                dSE2xdR = -sind(RW) * S_E2m(1) + cosd(RW) * S_E2m(3);
                dSE2zdR = -cosd(RW) * S_E2m(1) - sind(RW) * S_E2m(3);
                dSF2xdR = -sind(RW) * S_F2m(1) + cosd(RW) * S_F2m(3);
                dSF2zdR = -cosd(RW) * S_F2m(1) - sind(RW) * S_F2m(3);
            
                df1 = dSE2zdR * (CoR_Lz + S_E2(2) - E1(2));
                df2 = dSF2zdR * (CoR_Lz + S_F2(2) - F1(2));
                df = -dSE2xdR - 0.5 * df1 / f2 + dSF2xdR + 0.5 * df2 / f4;
            
                RW = RW - f/df;

                %Ensuring joint angle is within one rotation
                if RW < 0 | RW > 360
                    RW = wrapTo360(real(RW));
                else
                    RW = real(RW);
                end
            
                %Ensuring only one possible solution
                if RW < 180 && RW > 90
                    RW = RW + 180;
                elseif RW < 270 && RW > 180
                    RW = RW - 180;
                end
            
                Ay = A_rot(0,RW,0);
                S_E2 = Ay*S_E2m;
                S_F2 = Ay*S_F2m;
            
                f1_new = E1(1) - S_E2(1);
                f2_new = sqrt(L_short^2 - (CoR_Ly+S_E2(2)-E1(2))^2 - (CoR_Lz+S_E2(3)-E1(3))^2);
                f3_new = F1(1) - S_F2(1);
                f4_new = sqrt(L_short^2 - (CoR_Ly+S_F2(2)-F1(2))^2 - (CoR_Lz+S_F2(3)-F1(3))^2);
                f_new = f1_new - f2_new - f3_new + f4_new;
            
                errorW = abs(f_new);
                CountW = CountW + 1;
            end
            %Finding vectors from CoR to joints in global coordinate system
            S_C2 = A_rot(0,RE,0)*S_C2m;
            S_D2 = A_rot(0,RE,0)*S_D2m;
            S_E2 = A_rot(0,RW,0)*S_E2m;
            S_F2 = A_rot(0,RW,0)*S_F2m;
            
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
            
            %Test for length of linkages
            L_C = norm(C2-C1);
            L_D = norm(D2-D1);
            L_E = norm(E2-E1);
            L_F = norm(F2-F1);
            Test_C = SolutionCheck(L_C, L_short, 0.001);
            Test_D = SolutionCheck(L_D, L_short, 0.001);
            Test_E = SolutionCheck(L_E, L_short, 0.001);
            Test_F = SolutionCheck(L_F, L_short, 0.001);

            %Finding distances between vectors
            v_A = A2 - A1;
            v_B = B2 - B1;
            v_C = C2 - C1;
            v_D = D2 - D1;
            v_E = E2 - E1;
            v_F = F2 - F1;

            d_AF = D_lineline(A1,v_A,F1,v_F);
            d_BC = D_lineline(B1,v_B,C1,v_C);
            d_CD = D_lineline(C1,v_C,D1,v_D);
            d_EF = D_lineline(E1,v_E,F1,v_F);

            %Checking for irrelevant results
            if ~isreal(pz)
                break
            elseif d_AF < d | d_BC < d | d_CD < d | d_EF < d
                break
            end

            %Plotting if NR solver is within error
            if breaksignal == boolean(0) && Test_C && Test_D && Test_D && Test_E && Test_F
                Resultsx_neg(i,j) = -px;
                JointCoordx_negE(i,j) = RE;
                JointCoordx_negW(i,j) = RW;
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

figure(2)
surf(AxesY,AxesZ,Resultsx_neg)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('rotation [deg]')
title('Max rotation negative X')


figure(3)
subplot(3,1,1)
surf(AxesY,AxesZ,Resultsx_pos)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('rotation [deg]')
title('Max rotation positive X')
sgtitle('Positive X rotation')

subplot(3,1,2)
surf(AxesY,AxesZ,JointCoordx_posE)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Joint coordinate C [deg]')
title('Maximum Joint Coordinate E')

subplot(3,1,3)
surf(AxesY,AxesZ,JointCoordx_posW)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Joint coordinate W [deg]')
title('Maximum Joint Coordinate W')

figure(4)
subplot(3,1,1)
surf(AxesY,AxesZ,Resultsx_neg)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('rotation [deg]')
title('Max rotation negative X')
sgtitle('Negative X rotation')

subplot(3,1,2)
surf(AxesY,AxesZ,JointCoordx_negE)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Joint coordinate E [deg]')
title('Maximum Joint Coordinate E')

subplot(3,1,3)
surf(AxesY,AxesZ,JointCoordx_negW)
grid on
colorbar
xlabel('Y-axis [mm]')
ylabel('Z-axis [mm]')
zlabel('Joint coordinate W [deg]')
title('Maximum Joint Coordinate W')

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
plot3(E1(1),E1(2),E1(3),'.' ,E2(1), E2(2), E2(3),'.','MarkerSize',15, 'color','k')
text(E1(1)-0.05,E1(2)+0.1,E1(3)+0.05,'E1', 'FontSize',8)
text(E2(1)+0.05,E2(2)+0.05,E2(3)+0.05,'E2', 'FontSize',8)
plot3(F1(1),F1(2),F1(3),'.' ,F2(1), F2(2), F2(3), '.','MarkerSize',15, 'color','k')
text(F1(1)-0.05,F1(2)+0.1,F1(3)+0.05,'F1', 'FontSize',8)
text(F2(1)+0.05,F2(2)+0.05,F2(3)+0.05,'F2', 'FontSize',8)
plot3(TCP(1),TCP(2),TCP(3), '.','MarkerSize',15, 'color','k')
text(TCP(1)-0.05,TCP(2)+0.15,TCP(3)-0.05,'TCP', 'FontSize',8)
plot3(CoR_L(1),CoR_L(2),CoR_L(3),'.' ,CoR_R(1), CoR_R(2), CoR_R(3), '.','MarkerSize',15, 'color','k')
text(CoR_L(1)-0.05,CoR_L(2)+0.15,CoR_L(3)-0.05,'CoR_L', 'FontSize',8)
text(CoR_R(1)-0.05,CoR_R(2)+0.15,CoR_R(3)-0.05,'CoR_R', 'FontSize',8)

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

function y = SolutionCheck(Actual, Check, Range)
    %This function wil find if the solution the NR solver has found is
    %within a range. It will return boolean 1 if it is and boolean 0 if not
    if Actual < Check - Range | Actual > Check + Range
        y = boolean(0);
    else
        y = boolean(1);
    end
end