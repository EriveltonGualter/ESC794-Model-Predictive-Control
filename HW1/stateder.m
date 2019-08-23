%state derivatives
%4 -dof robot (no gravity effects)

function zdot=stateder(z,u)
%z is the stacked vector of all joint positions (z1) and all joint
%velocities (z2)

%assign q, qdot from z

z1=z(1:4);
z2=z(5:8);

q1=z1(1);q2=z1(2);q3=z1(3);q4=z1(4);
q1dot=z2(1);q2dot=z2(2);q3dot=z2(3);q4dot=z2(4);

%Imported symbolic code

M=[ 0.07214*sin(2.0*q2) + 0.01124*sin(2.0*q4) + 0.3504*cos(q4) + 0.1056*sin(q4) - 0.9447*cos(q2)^2 + 0.0504*cos(q4)^2 - 0.359*cos(q2)^2*cos(q4) - 0.07698*cos(q2)^2*sin(q4) - 0.1008*cos(q2)^2*cos(q4)^2 - 0.04496*cos(q2)*cos(q4)^2*sin(q2) - 0.04496*cos(q2)^2*cos(q4)*sin(q4) - 0.07698*cos(q2)*cos(q4)*sin(q2) + 0.359*cos(q2)*sin(q2)*sin(q4) + 0.1008*cos(q2)*cos(q4)*sin(q2)*sin(q4) + 1.143, 5.199e-5*cos(q2)*cos(q4) - 0.0005887*sin(q2) - 0.003375*cos(q2) + 4.624e-6*cos(q2)*sin(q4) + 4.624e-6*cos(q4)*sin(q2) - 5.199e-5*sin(q2)*sin(q4), 0.06538*cos(q2) + 0.07231*sin(q2) - 0.008644*cos(q2)*cos(q4) + 0.02867*cos(q2)*sin(q4) - 0.03849*cos(q4)*sin(q2) + 0.1795*sin(q2)*sin(q4) - 0.0504*cos(q2)*cos(q4)^2 - 0.02248*cos(q4)^2*sin(q2) - 0.02248*cos(q2)*cos(q4)*sin(q4) + 0.0504*cos(q4)*sin(q2)*sin(q4), 5.199e-5*cos(q2)*cos(q4) + 4.624e-6*cos(q2)*sin(q4) + 4.624e-6*cos(q4)*sin(q2) - 5.199e-5*sin(q2)*sin(q4);
                                                                                                                                                                                                                                                5.199e-5*cos(q2)*cos(q4) - 0.0005887*sin(q2) - 0.003375*cos(q2) + 4.624e-6*cos(q2)*sin(q4) + 4.624e-6*cos(q4)*sin(q2) - 5.199e-5*sin(q2)*sin(q4),                                                                                                          0.3418*cos(q4) + 0.1343*sin(q4) + 1.077,                                                                                                                                                                                                                     5.199e-5*cos(q4) + 4.624e-6*sin(q4) + 0.0003001,                                                                0.1709*cos(q4) + 0.06716*sin(q4) + 0.06095;
                                                                                                                             0.06538*cos(q2) + 0.07231*sin(q2) - 0.008644*cos(q2)*cos(q4) + 0.02867*cos(q2)*sin(q4) - 0.03849*cos(q4)*sin(q2) + 0.1795*sin(q2)*sin(q4) - 0.0504*cos(q2)*cos(q4)^2 - 0.02248*cos(q4)^2*sin(q2) - 0.02248*cos(q2)*cos(q4)*sin(q4) + 0.0504*cos(q4)*sin(q2)*sin(q4),                                                                                                  5.199e-5*cos(q4) + 4.624e-6*sin(q4) + 0.0003001,                                                                                                                                                                           0.02867*sin(q4) - 0.008644*cos(q4) - 0.02248*cos(q4)*sin(q4) - 0.0504*cos(q4)^2 + 0.06538,                                                                       5.199e-5*cos(q4) + 4.624e-6*sin(q4);
                                                                                                                                                                                                                                                                                       5.199e-5*cos(q2)*cos(q4) + 4.624e-6*cos(q2)*sin(q4) + 4.624e-6*cos(q4)*sin(q2) - 5.199e-5*sin(q2)*sin(q4),                                                                                                       0.1709*cos(q4) + 0.06716*sin(q4) + 0.06095,                                                                                                                                                                                                                                 5.199e-5*cos(q4) + 4.624e-6*sin(q4),                                                                                                   0.06095];
 
C=[                                                    0.1795*q2dot*sin(2.0*q2 + q4) - 0.01924*q4dot*cos(2.0*q2 + q4) - 0.03849*q2dot*cos(2.0*q2 + q4) + 0.08976*q4dot*sin(2.0*q2 + q4) + 0.0609*q2dot*cos(2.0*q2) + 0.4976*q2dot*sin(2.0*q2) + 0.03358*q4dot*cos(q4) - 0.01124*q2dot*cos(2.0*q2 + 2.0*q4) - 0.01124*q4dot*cos(2.0*q2 + 2.0*q4) - 0.08544*q4dot*sin(q4) + 0.0252*q2dot*sin(2.0*q2 + 2.0*q4) + 0.0252*q4dot*sin(2.0*q2 + 2.0*q4), 0.1795*q1dot*sin(2.0*q2 + q4) - 0.01679*q3dot*cos(q2 - 1.0*q4) - 0.00562*q3dot*cos(q2 + 2.0*q4) - 0.03849*q1dot*cos(2.0*q2 + q4) - 0.04272*q3dot*sin(q2 - 1.0*q4) + 0.0126*q3dot*sin(q2 + 2.0*q4) + 0.0609*q1dot*cos(2.0*q2) + 0.4976*q1dot*sin(2.0*q2) + 4.624e-6*q2dot*cos(q2 + q4) - 0.002455*q3dot*cos(q2 + q4) + 4.624e-6*q4dot*cos(q2 + q4) - 5.199e-5*q2dot*sin(q2 + q4) + 0.04704*q3dot*sin(q2 + q4) - 5.199e-5*q4dot*sin(q2 + q4) - 0.0005887*q2dot*cos(q2) + 0.03053*q3dot*cos(q2) - 0.01124*q1dot*cos(2.0*q2 + 2.0*q4) + 0.003375*q2dot*sin(q2) - 0.02009*q3dot*sin(q2) + 0.0252*q1dot*sin(2.0*q2 + 2.0*q4), 0.01679*q4dot*cos(q2 - 1.0*q4) - 0.00562*q2dot*cos(q2 + 2.0*q4) - 0.01679*q2dot*cos(q2 - 1.0*q4) - 0.01124*q4dot*cos(q2 + 2.0*q4) - 0.04272*q2dot*sin(q2 - 1.0*q4) + 0.0126*q2dot*sin(q2 + 2.0*q4) + 0.04272*q4dot*sin(q2 - 1.0*q4) + 0.0252*q4dot*sin(q2 + 2.0*q4) - 0.002455*q2dot*cos(q2 + q4) - 0.002455*q4dot*cos(q2 + q4) + 0.04704*q2dot*sin(q2 + q4) + 0.04704*q4dot*sin(q2 + q4) + 0.03053*q2dot*cos(q2) - 0.02009*q2dot*sin(q2), 0.01679*q3dot*cos(q2 - 1.0*q4) - 0.01924*q1dot*cos(2.0*q2 + q4) - 0.01124*q3dot*cos(q2 + 2.0*q4) + 0.08976*q1dot*sin(2.0*q2 + q4) + 0.04272*q3dot*sin(q2 - 1.0*q4) + 0.0252*q3dot*sin(q2 + 2.0*q4) + 4.624e-6*q2dot*cos(q2 + q4) - 0.002455*q3dot*cos(q2 + q4) + 4.624e-6*q4dot*cos(q2 + q4) - 5.199e-5*q2dot*sin(q2 + q4) + 0.04704*q3dot*sin(q2 + q4) - 5.199e-5*q4dot*sin(q2 + q4) + 0.03358*q1dot*cos(q4) - 0.01124*q1dot*cos(2.0*q2 + 2.0*q4) - 0.08544*q1dot*sin(q4) + 0.0252*q1dot*sin(2.0*q2 + 2.0*q4);
 0.03849*q1dot*cos(2.0*q2 + q4) + 0.01679*q3dot*cos(q2 - 1.0*q4) + 0.00562*q3dot*cos(q2 + 2.0*q4) - 0.1795*q1dot*sin(2.0*q2 + q4) + 0.04272*q3dot*sin(q2 - 1.0*q4) - 0.0126*q3dot*sin(q2 + 2.0*q4) - 0.0609*q1dot*cos(2.0*q2) - 0.4976*q1dot*sin(2.0*q2) + 0.002455*q3dot*cos(q2 + q4) - 0.04704*q3dot*sin(q2 + q4) - 0.03053*q3dot*cos(q2) + 0.01124*q1dot*cos(2.0*q2 + 2.0*q4) + 0.02009*q3dot*sin(q2) - 0.0252*q1dot*sin(2.0*q2 + 2.0*q4),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               q4dot*(0.06716*cos(q4) - 0.1709*sin(q4)),                                                                                                                                                                                         q1dot*(0.01679*cos(q2 - 1.0*q4) + 0.00562*cos(q2 + 2.0*q4) + 0.04272*sin(q2 - 1.0*q4) - 0.0126*sin(q2 + 2.0*q4) + 0.002455*cos(q2 + q4) - 0.04704*sin(q2 + q4) - 0.03053*cos(q2) + 0.02009*sin(q2)) + q4dot*(2.312e-6*cos(q4) - 2.599e-5*sin(q4)),                                                                                                                                                                                                                                                                                                                                                                              q2dot*(0.06716*cos(q4) - 0.1709*sin(q4)) + q4dot*(0.06716*cos(q4) - 0.1709*sin(q4)) + q3dot*(2.312e-6*cos(q4) - 2.599e-5*sin(q4));
   0.01679*q4dot*cos(q2 - 1.0*q4) - 0.00562*q2dot*cos(q2 + 2.0*q4) - 0.01679*q2dot*cos(q2 - 1.0*q4) - 0.01124*q4dot*cos(q2 + 2.0*q4) - 0.04272*q2dot*sin(q2 - 1.0*q4) + 0.0126*q2dot*sin(q2 + 2.0*q4) + 0.04272*q4dot*sin(q2 - 1.0*q4) + 0.0252*q4dot*sin(q2 + 2.0*q4) - 0.002455*q2dot*cos(q2 + q4) - 0.002455*q4dot*cos(q2 + q4) + 0.04704*q2dot*sin(q2 + q4) + 0.04704*q4dot*sin(q2 + q4) + 0.03053*q2dot*cos(q2) - 0.02009*q2dot*sin(q2),                                                                                                                                                                                                                                                                                                                                                                  q4dot*(2.312e-6*cos(q4) - 2.599e-5*sin(q4)) - 1.0*q1dot*(0.01679*cos(q2 - 1.0*q4) + 0.00562*cos(q2 + 2.0*q4) + 0.04272*sin(q2 - 1.0*q4) - 0.0126*sin(q2 + 2.0*q4) + 0.002455*cos(q2 + q4) - 0.04704*sin(q2 + q4) - 0.03053*cos(q2) + 0.02009*sin(q2)),                                                                                                                                                                                                                                                                                                                                         q4dot*(0.01433*cos(q4) + 0.004322*sin(q4) + 0.0504*cos(q4)*sin(q4) - 0.02248*cos(q4)^2 + 0.01124),                                                                                                        0.01679*q1dot*cos(q2 - 1.0*q4) - 0.01124*q1dot*cos(q2 + 2.0*q4) + 0.04272*q1dot*sin(q2 - 1.0*q4) + 0.0252*q1dot*sin(q2 + 2.0*q4) - 0.01124*q3dot*cos(2.0*q4) + 0.0252*q3dot*sin(2.0*q4) - 0.002455*q1dot*cos(q2 + q4) + 0.04704*q1dot*sin(q2 + q4) + 2.312e-6*q2dot*cos(q4) + 0.01433*q3dot*cos(q4) + 4.624e-6*q4dot*cos(q4) - 2.599e-5*q2dot*sin(q4) + 0.004322*q3dot*sin(q4) - 5.199e-5*q4dot*sin(q4);
                                                      0.01924*q1dot*cos(2.0*q2 + q4) - 0.01679*q3dot*cos(q2 - 1.0*q4) + 0.01124*q3dot*cos(q2 + 2.0*q4) - 0.08976*q1dot*sin(2.0*q2 + q4) - 0.04272*q3dot*sin(q2 - 1.0*q4) - 0.0252*q3dot*sin(q2 + 2.0*q4) + 0.002455*q3dot*cos(q2 + q4) - 0.04704*q3dot*sin(q2 + q4) - 0.03358*q1dot*cos(q4) + 0.01124*q1dot*cos(2.0*q2 + 2.0*q4) + 0.08544*q1dot*sin(q4) - 0.0252*q1dot*sin(2.0*q2 + 2.0*q4),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       - 1.0*q2dot*(0.06716*cos(q4) - 0.1709*sin(q4)) - 1.0*q3dot*(2.312e-6*cos(q4) - 2.599e-5*sin(q4)),                                                                                     0.01124*q1dot*cos(q2 + 2.0*q4) - 0.01679*q1dot*cos(q2 - 1.0*q4) - 0.04272*q1dot*sin(q2 - 1.0*q4) - 0.0252*q1dot*sin(q2 + 2.0*q4) + 0.01124*q3dot*cos(2.0*q4) - 0.0252*q3dot*sin(2.0*q4) + 0.002455*q1dot*cos(q2 + q4) - 0.04704*q1dot*sin(q2 + q4) - 2.312e-6*q2dot*cos(q4) - 0.01433*q3dot*cos(q4) + 2.599e-5*q2dot*sin(q4) - 0.004322*q3dot*sin(q4),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              0];
                                                  
                                                  
 %State derivatives
 z1dot=z2;
 z2dot=inv(M)*(-C*z2+u);
 zdot=[z1dot;z2dot];
 