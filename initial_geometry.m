%{
initial_geometry.m

Created: 11/15/2020
Last Updated: 05/09/2023

written by Andrew Schmidt
Cardiomyocyte Laboratory
Henry Samueli School of Engineering
University of California, Irvine, Irvine, CA 92697

The purpose of this code is to generate the initial geometries used in the SarcContractionCycle

Input: 
        
Output: .mat file of all the arrays describing the initial node placement (geometry) 
        for the half-sarcomere contraction simulation
%}

function initial_geometry()
%% Define Geometrical Constants
%Define number of nodes and rest lengths
AMpairs = [1, 1; 9, 6; 18, 12; 24, 16; 32, 20; 40, 25; 48, 30; 64, 40; 80, 50; 320, 200];
t_ends =    [2, 2.5,   4.5,     6,     7.5,    8.5,     11,     14,     18, 300];
caseno = 9;
t_end = t_ends(caseno);
a_nodes = AMpairs(caseno, 1); %number of actin nodes
m_nodes = AMpairs(caseno, 2); %number of myosin nodes (aka XBs)
int_nodes = a_nodes + m_nodes; %number of interior nodes (not Z- or M-line)
end_nodes = 2; %Z-line and M-line nodes
tot_nodes = int_nodes + end_nodes;
m0 = 14.3e-9; %42.9e-9; %rest length of myosin spring (distance between myosin nodes)14.3 nm in m (Chase, 2004)
a0 = 5.5e-9; %5.5e-9; %36.9e-9; %m0/2; %12.3e-9; %12.3e-9; %rest length of actin spring (distance between actin nodes) 12.3 nm in m (Chase, 2004)
xb0 = 10e-9; %rest length of XB spring 12 nm extended in m <-- probably depends on spring constant?**********
%xhs = 43e-9; %rest length of half sarcomere 43 nm in m*********************** redefined for mechanical equilibrium later
xf0 = 100e-9; %rest length of external force spring (representing pdms)
%xf = xhs + xf0; %location of external force spring anchor
%titin0 = xhs; %in m when titin connects Z-line to M-line
%titin0 = xhs - m_nodes*m_0 - az; %in m when titin connects Z-line to m1
d_ps = 7e-9; %7 nm in m (Chase, 2004)
%titin0 = currently defined in test cases


%Define time parameters, this will not give the total time that the simulation will undergo contraction,
%it will only give the total number of time steps the simulation will undergo.
%The total time will depend on the time step values calculated in the main code.
%The time step values are variable depending on the current state of the system.
dt = 1e-3; %Maximum timestep size, in s
tot_timesteps = round(t_end/dt)+1;


N_locations = zeros(tot_timesteps, tot_nodes); %rows = time, columns are locations of all nodes at each time (incl. Z- and M-line)
%% Test Cases
%*************************************************************************
%Initialize node locations
N_locations(1, 1:(a_nodes+1)) = a0.*[0:1:a_nodes]; %a_positions
testcase = 0;
% tc_suffix = 0.1*a_nodes + 0.01*m_nodes; %Works as long as a_nodes and m_nodes both < 10
%{
%Normal(?) node locations
% N_locations(1, (a_nodes+2:tot_nodes)) = xhs-m0.*[m_nodes:-1:0]; %m_positions
N_locations(1, (a_nodes+2:tot_nodes)) = N_locations(1, 2) + xb0.*[1:1:m_nodes+1]; %m_positions
xf = N_locations(1, tot_nodes)+xf0;
titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
%}


%TEST CASE 1: Example where there displacement x - xb0 = 0 between actin and myosin nodes at time = 0; m0=a0
%{
testcase = 1;
testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];

m0 = a0;
N_locations(1, ((a_nodes+2):(tot_nodes-1)))= N_locations(1, 2:a_nodes+1)+ xb0; %a0.*[1:1:m_nodes]+xb0;
N_locations(1, tot_nodes) = N_locations(1, tot_nodes-1) + m0;
xf = N_locations(1, tot_nodes)+xf0;
%titin0 = N_locations(1, tot_nodes) - N_locations(1, 1); %titin connects Z-line to M-line
titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
%}

%TEST CASE 2: Only myosin 5 perfectly aligned, myosin 4 is m0 away from myosin 5; m0=a0
%{
testcase = 2;
testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];

m0 = a0;
N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ xb0; %myosin 5 location
N_locations(1, (a_nodes+2)) = N_locations(1, (tot_nodes-1)) - m0; %myosin 4 location
N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0;
xf = N_locations(1, tot_nodes)+xf0;
%titin0 = N_locations(1, tot_nodes) - N_locations(1, 1); %titin connects Z-line to M-line
titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
%}

%% TEST CASE 3: Perfectly aligned; m0 ~= a0 
%{
testcase = 3;
testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
%%{
if a_nodes == 1 && m_nodes == 1
    N_locations(1, (a_nodes+2)) = N_locations(1, a_nodes+1)+ xb0; %myosin 1 location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
    xf = N_locations(1, tot_nodes) + xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
end
%}
%{
if a_nodes == 1 && m_nodes == 1 %offset to match 6.62_3pullbias
    N_locations(1, (a_nodes+2)) = N_locations(1, a_nodes+1) + xb0 + m0/4 - 0.5e-9; %myosin 1 location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
    xf = N_locations(1, tot_nodes) + xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    testcase = [testcase, '_pullbias'];
end
%}  
    
%}

%TEST CASE 3 VERSION 2: varied distances from perfectly aligned
%{
for offset = -5e-9:1e-9:5e-9 %in m
    testcase = 3;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes), '_offset=', num2str(offset*1e9), 'nm'];

    if a_nodes == 1 && m_nodes == 1
        N_locations(1, (a_nodes+2)) = N_locations(1, a_nodes+1) + xb0 + offset; %myosin 1 location
        N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
        xf = N_locations(1, tot_nodes) + xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
end %<- move to after "save(filename)" when running this version of TEST CASE 3
%}

%% TEST CASE 4: Neither is perfectly aligned; a0 = 12.3 nm, m0 = 14.3 nm
%{
testcase = 4;
testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];

if m_nodes == 2 && a_nodes == 2
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1) + m0; %myosin 5 location
    N_locations(1, (a_nodes+2)) = N_locations(1, (tot_nodes-1)) - m0; %myosin 4 location
    
    %Alternative: myosin 5 is m0 in front of first actin, not second actin binding site
%     N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes)+ m0; %myosin 5 location
%     N_locations(1, (a_nodes+2)) = N_locations(1, (tot_nodes-1)) - m0; %myosin 4 location
    
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    
elseif a_nodes == 1 && m_nodes == 1
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ m0; %myosin location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    
elseif a_nodes > m_nodes
    % Final myosin is m0 in front of FINAL actin node, rest of myosin use final myosin as reference
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ m0; %myosin closest to M-line location
    N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
    N_locations(1, tot_nodes) = N_locations(1, tot_nodes-1) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    
elseif a_nodes == 4 && m_nodes == 4
    % Final myosin is 2*m0 in front of SECOND TO LAST actin node
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes)+ 2*m0; %myosin closest to M-line location
    N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
    N_locations(1, tot_nodes) = N_locations(1, tot_nodes-1) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    
else %Any number of m_nodes and a_nodes as long as m_nodes >= a_nodes
    % Final myosin is m0 in front of FINAL actin node, rest of myosin use final myosin as reference
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ m0; %myosin closest to M-line location
    N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
    N_locations(1, tot_nodes) = N_locations(1, tot_nodes-1) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
end
%}

%% TEST CASE 5: Both are perfectly aligned; a0 = 12.3 nm, m0 = a0 nm %PICKUP
%{
testcase = 5;
%Uncomment only one of the following 5.22_#

%TC 5.22_1 Both myosin are perfectly aligned with their actin sites at the start
% if m_nodes == 2 && a_nodes == 2 
%     testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes), '_1'];
% %     m0 = a0;
%     N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ xb0; %myosin 2 location
%     N_locations(1, (a_nodes+2)) = N_locations(1, (a_nodes)) + xb0; %myosin 1 location
%     N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %m-line location
% %     N_locations(1, tot_nodes) = N_locations(1, a_nodes+2) + m0*m_nodes;
%     xf = N_locations(1, tot_nodes) + xf0;
%     titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
%     xhs = N_locations(1, tot_nodes) - N_locations(1, 1);
% end

%TC 5.22_2 m1 perfectly aligned at t = 0, when m1 enters state 3 (power stroke), m2 is perfectly aligned to bind to actin site
testcase = 5;
if m_nodes == 2 && a_nodes == 2
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes), '_2'];
%     m0 = a0;
    N_locations(1, (a_nodes+2)) = N_locations(1, (a_nodes)) + xb0; %myosin 1 location
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+2) + m0; %myosin 2 location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %m-line location
    xf = N_locations(1, tot_nodes) + xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1);
end
%}

%% TEST CASE 6: Neither is perfectly aligned; a0 = 12.3 nm, m0 = 14.3 nm, a_nodes > m_nodes, a_nodes > 2
%{
    testcase = 6;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)]; % Final myosin is m0 in front of SECOND actin, rest of myosin use final myosin as reference
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ m0; %myosin closest to M-line location
    
    %TEST CASE 6.52 --> 5 actin a0 = 5.5 nm 2 myosin
    % Final myosin is m0 in front of THIRD actin, rest of myosin use final myosin as reference
    if a_nodes == 5 && m_nodes == 2 && a0 == 5.5e-9
        N_locations(1, (tot_nodes-1)) = N_locations(1, 4)+ m0; %myosin closest to M-line location
    end
    
    %TEST CASE 6.62 --> 6 actin a0 = 5.5 nm 2 myosin
    % Final myosin is m0 - 1.55e-9 in front of FOURTH actin, rest of myosin use final myosin as reference
    if a_nodes == 6 && m_nodes == 2 && a0 == 5.5e-9
        N_locations(1, (tot_nodes-1)) = N_locations(1, 5) + m0 - 1.55e-9; %myosin closest to M-line location centered between 4th and 5th actin
    end
  % Uncomment only one of the following since outside the geometry they are the same conditions set
    % TC 6.62_1 Final myosin is xb0 in front of FOURTH actin, rest of myosin use final myosin as reference
%     if a_nodes == 6 && m_nodes == 2 && a0 == m0/2
%         N_locations(1, (tot_nodes-1)) = N_locations(1, 5) + xb0; %myosin closest to M-line location centered between 4th and 5th actin
%     end %changed to xb0 instead of m0 in front of fourth actin on 10/2/2021
    % TC 6.62_2 Final myosin is m0/4 + xb0 in front of FOURTH actin, rest of myosin use final myosin as reference - end of XB are equidistant from two nearest actin
%     if a_nodes == 6 && m_nodes == 2 && a0 == m0/2
%         N_locations(1, (tot_nodes-1)) = N_locations(1, 5) + (m0/4 + xb0); %myosin closest to M-line location centered between 4th and 5th actin
%     end
    % TC 6.62_3_pullbias Final myosin is m0 in front of FOURTH actin, rest of myosin use final myosin as reference
    if a_nodes == 6 && m_nodes == 2 && a0 == m0/2
        N_locations(1, (tot_nodes-1)) = N_locations(1, 5) + (m0/4 - 0.5e-9 + xb0); %myosin closest to M-line location off-center between 4th and 5th actin
        testcase = [testcase, '_3pullbias'];
    end
    
    % TC 6.62_3_phys Final myosin is a barezone of 3.35 nm from m-line and total length of half sarcomere is 51.6 nm
%     N_locations(1, tot_nodes) = 51.6e-9; %<- have to replace M-line location at end of section with this line too
%     if a_nodes == 6 && m_nodes == 2 && a0 == m0/2
%         N_locations(1, (tot_nodes-1)) = 51.6e-9 - 3.35e-9; %last myosin is a 3.35 nm barezone away from m-line
%         testcase = [testcase, '_3phys'];
%     end
    
    %TEST CASE 6.64 --> 6 actin, a0 = 12.3 nm 4 myosin
    if a_nodes == 6 && m_nodes == 4
        %Final m_nodes is m0 in front of SECOND TO LAST actin node 
        N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes)+ m0; %myosin closest to M-line location
    end
    
    %TEST CASE 6.46 --> 4 actin, 6 myosin
    if a_nodes == 4 && m_nodes == 6
        %Final m_nodes is 3*m0 in front of LAST actin node 
        N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+1)+ 3*m0; %myosin closest to M-line location
    end
    
    %TEST CASE 6.42 --> 4 actin, 2 myosin
    if a_nodes == 4 && m_nodes == 2
        %Final m_nodes is m0 in front of SECOND TO LAST actin node 
        N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes)+ m0; %myosin closest to M-line location
    end
    
    
    
    %Rest of test 6.__ cases
    N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; %m_nodes locations
    N_locations(1, tot_nodes) = N_locations(1, tot_nodes-1) + m0;
    xf = N_locations(1, tot_nodes)+xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node    
%}

%% TEST CASE 7: Many (any number) actin nodes and 2 myosin nodes; a0 = 12.3 m0 = 14.3;
%{
% m0 = a0;
testcase = 7;
testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
if m_nodes == 2
    N_locations(1, (a_nodes+2)) = N_locations(1, a_nodes+1) + xb0; %myosin 4 location
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes+2)+ m0; %myosin 5 location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
    xf = N_locations(1, tot_nodes) + xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
else %if a0 == m0/2
%     N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes) + (m0/4 - 0.5e-9 + xb0); %[when a0 = m0/2] myosin closest to M-line location is off-center between last and second to last actin
    N_locations(1, (tot_nodes-1)) = N_locations(1, a_nodes) + (a0/2 - 0.5e-9 + xb0); %[when a0 = 12.3 nm] myosin closest to M-line location is off-center between last and second to last actin
    N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0];
%     N_locations(1, (a_nodes+2)) = N_locations(1, tot_nodes-1) - m0; %myosin 1 location
    N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
    xf = N_locations(1, tot_nodes) + xf0;
    titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
    xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
end
%}


%% TEST CASE 8: Same conditions as 6, but with many actin and many myosin: No myosin is perfectly aligned; a0 = 12.3 nm, m0 = 14.3 nm, a_nodes > m_nodes, a_nodes >> 2
%%{
    testcase = 8;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
    %TEST CASE 8.32 --> 3 actin, 2 myosin
    if a_nodes == 3 && m_nodes == 2
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 6.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.96 --> 9 actin, 6 myosin 
    if a_nodes == 9 && m_nodes == 6
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 18.14e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.1812 --> 18 actin, 12 myosin
    if a_nodes == 18 && m_nodes == 12
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 36.39e-9;  %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.2416 --> 24 actin, 16 myosin, Barezone of 24.37 nm between last myosin node and M-line
    if a_nodes == 24 && m_nodes == 16
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 48.37e-9; %59.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %24.37e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    if a_nodes == 24 && m_nodes == 16 && a0 == 5.5e-9 %Geeves and Holmes Actin monomer spacing on single actin filament
        %Update values
        a_nodes = 54; %a_nodes = 24*12.3nm / 5.5 nm
        int_nodes = a_nodes + m_nodes; %number of interior nodes (not Z- or M-line)
        tot_nodes = int_nodes + end_nodes;
        N_locations = zeros(tot_timesteps, tot_nodes);
        N_locations(1, 1:(a_nodes+1)) = a0.*[0:1:a_nodes];
        testcase = 8;
        testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
        
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 48.37e-9; %59.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %24.37e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    if a_nodes == 24 && m_nodes == 16 && a0 == 2.7e-9 %Geeves and Holmes Actin monomer spacing on double filament
        %Update values
        a_nodes = 110; %a_nodes = 24*12.3nm / 2.7 nm
        int_nodes = a_nodes + m_nodes; %number of interior nodes (not Z- or M-line)
        tot_nodes = int_nodes + end_nodes;
        N_locations = zeros(tot_timesteps, tot_nodes);
        N_locations(1, 1:(a_nodes+1)) = a0.*[0:1:a_nodes];
        testcase = 8;
        testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
        
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 48.37e-9; %59.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %24.37e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.3220 --> 32 actin, 20 myosin
    if a_nodes == 32 && m_nodes == 20
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 60.65e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.4025 --> 1/2 HALF SARCOMERE 40 actin, 25 myosin, Barezone of 38.5 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 40 && m_nodes == 25
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 75.8e-9; %100e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %38.5e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
        
    %TEST CASE 8.4830 --> 48 actin, 30 myosin
    if a_nodes == 48 && m_nodes == 30
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 90.98e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.6440 --> 1/2 HALF SARCOMERE 60 actin, 40 myosin, Barezone of 121.3 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 64 && m_nodes == 40
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 121.3e-9; %100e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %38.5e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.8050 --> FULL HALF SARCOMERE 80 actin, 50 myosin, Barezone of 77 nm between last myosin node and M-line, 200 nm between last actin and M-line
    if a_nodes == 80 && m_nodes == 50
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 151.6e-9; %137.3e-9; <- what was being used prior to 2/4/2022 %200e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %61.6e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    if a_nodes == 80 && m_nodes == 50 && a0 == 5.5e-9 %Geeves and Holmes Actin monomer spacing on single actin filament
        %Update values
        a_nodes = ceil(a_nodes*12.3e-9/a0); %a_nodes = 24*12.3nm / 5.5 nm
        int_nodes = a_nodes + m_nodes; %number of interior nodes (not Z- or M-line)
        tot_nodes = int_nodes + end_nodes;
        clear N_locations
        N_locations = zeros(tot_timesteps, tot_nodes);
        N_locations(1, 1:(a_nodes+1)) = a0.*[0:1:a_nodes];
        testcase = 8;
        testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
        
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + m0; %59.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %24.37e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.101 --> Goal: both m1 and m2 have the same duty ratio across different conditions
    if a_nodes == 10 && m_nodes == 2 %INCOMPLETE
        a0 = 5.5e-9;
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 137.3e-9; %200e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %77e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.320200 --> Very long sarcomere, used to evaluate whether Tau trends persists at large myosin number values
    if a_nodes == 320 && m_nodes == 200
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 606.5e-9; %137.3e-9;
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end

%}

%% TEST CASE 9: Same conditions as 8, all sarcomeres the same length, the first actin and final myosin spring are very large to compensate for this; a0 = 12.3 nm, m0 = 14.3 nm, a_nodes > m_nodes, a_nodes >> 2
%{
    testcase = 9;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
    
    m_cond = 4; % 1 = 8.96, 2 = 8.2416, 3 = 8.4025, 4 = 8.6440
%     a_active_vec = [31, 39; 26, 50; 17, 56; 11, 74]; %active actin for 6 myosin, 16 myosin, 25 myosin, 40 myosin
    a_active_vec = [31, 26, 17, 11]; %start index of actin for 6 myosin, 16 myosin, 25 myosin, 40 myosin
    m_active_vec = [6, 16, 25, 40];
%     m_deact = [a_nodes+2+m_active_vec(m_cond):1:tot_nodes-1];
%     a_deact = [2:1:a_active_vec(m_cond, 1), a_active_vec(m_cond, 2):1:a_nodes+1];
    
    
    %initial actin nodes
%     N_locations(1, 1) = a0*a_active_vec(m_cond);
    N_locations(1, 2:(a_nodes+1)) = a0.*[a_active_vec(m_cond):1:a_active_vec(m_cond)+a_nodes-1]; %a_positions
    a0_long = a0*a_active_vec(m_cond);
    m0_long = (50-m_nodes+1)*m0;
    MLine = 1135.6e-9; %location of M-line in testcase 8.8050
    N_locations(1, tot_nodes) = MLine; %M-line location
    
    %TEST CASE 9.32 --> 3 actin, 2 myosin 
    if a_nodes == 3 && m_nodes == 2
%         actin_32 = [33, 35];
        N_locations(1, 2:(a_nodes+1)) = a0.*[33:1:35]; %a_positions
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; %rest of myosin
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end    
    
    %TEST CASE 9.96 --> 9 actin, 6 myosin 
    if a_nodes == 9 && m_nodes == 6
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    
    %TEST CASE 9.2416 --> 24 actin, 16 myosin, Barezone of 24.37 nm between last myosin node and M-line
    if a_nodes == 24 && m_nodes == 16
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
    %TEST CASE 9.4025 --> 1/2 HALF SARCOMERE 40 actin, 25 myosin, Barezone of 38.5 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 40 && m_nodes == 25
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
        
    %TEST CASE 9.6440 
    if a_nodes == 64 && m_nodes == 40
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
    %TEST CASE 9.8050 --> FULL HALF SARCOMERE 80 actin, 50 myosin, Barezone of 77 nm between last myosin node and M-line, 200 nm between last actin and M-line
    if a_nodes == 80 && m_nodes == 50
       N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
%}


%% TEST CASE 10: Same conditions as 8, Implement Bare Zone; a0 = 12.3 nm, m0 = 14.3 nm, a_nodes > m_nodes, a_nodes >> 2
%{
    testcase = 10;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
    
    m_cond = 4; % 1 = 8.96, 2 = 8.2416, 3 = 8.4025, 4 = 8.6440
%     a_active_vec = [31, 39; 26, 50; 17, 56; 11, 74]; %active actin for 6 myosin, 16 myosin, 25 myosin, 40 myosin
    a_active_vec = [31, 26, 17, 11]; %start index of actin for 6 myosin, 16 myosin, 25 myosin, 40 myosin
    m_active_vec = [6, 16, 25, 40];
%     m_deact = [a_nodes+2+m_active_vec(m_cond):1:tot_nodes-1];
%     a_deact = [2:1:a_active_vec(m_cond, 1), a_active_vec(m_cond, 2):1:a_nodes+1];
    
    
    %initial actin nodes
%     N_locations(1, 1) = a0*a_active_vec(m_cond);
    N_locations(1, 2:(a_nodes+1)) = a0.*[a_active_vec(m_cond):1:a_active_vec(m_cond)+a_nodes-1]; %a_positions
    a0_long = a0*a_active_vec(m_cond);
    m0_long = (50-m_nodes+1)*m0;
    MLine = 1135.6e-9; %location of M-line in testcase 8.8050
    N_locations(1, tot_nodes) = MLine; %M-line location
    
    %TEST CASE 9.32 --> 3 actin, 2 myosin 
    if a_nodes == 3 && m_nodes == 2
%         actin_32 = [33, 35];
        N_locations(1, 2:(a_nodes+1)) = a0.*[33:1:35]; %a_positions
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; %rest of myosin
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end    
    
    %TEST CASE 9.96 --> 9 actin, 6 myosin 
    if a_nodes == 9 && m_nodes == 6
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    
    %TEST CASE 9.2416 --> 24 actin, 16 myosin, Barezone of 24.37 nm between last myosin node and M-line
    if a_nodes == 24 && m_nodes == 16
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
    %TEST CASE 9.4025 --> 1/2 HALF SARCOMERE 40 actin, 25 myosin, Barezone of 38.5 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 40 && m_nodes == 25
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
        
    %TEST CASE 9.6440 
    if a_nodes == 64 && m_nodes == 40
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
    %TEST CASE 9.8050 --> FULL HALF SARCOMERE 80 actin, 50 myosin, Barezone of 77 nm between last myosin node and M-line, 200 nm between last actin and M-line
    if a_nodes == 80 && m_nodes == 50
       N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0_long; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
%}

%% TEST CASE 11: Same conditions as 8, Implement shift of some distance relative to the standard conditions in testcase 8 systems
%{

    testcase = 11;
    testcase = [num2str(testcase), '.', num2str(a_nodes), num2str(m_nodes)];
    
    shift = -0.5*m0;
    
    %TEST CASE 8.11 --> 1 actin, 1 myosin
    if a_nodes == 1 && m_nodes == 1
        N_locations(1, (a_nodes+2)) = N_locations(1, a_nodes+1)+ xb0 + shift; %myosin 1 location
        N_locations(1, tot_nodes) = N_locations(1, (tot_nodes-1)) + m0; %M-line location
        xf = N_locations(1, tot_nodes) + xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node
    end
    
    %TEST CASE 8.32 --> 3 actin, 2 myosin
    if a_nodes == 3 && m_nodes == 2
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 6.04e-9 + shift; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.96 --> 9 actin, 6 myosin 
    if a_nodes == 9 && m_nodes == 6
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 18.14e-9 + shift; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) - [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.1812 --> 18 actin, 12 myosin
    if a_nodes == 18 && m_nodes == 12
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 36.39e-9 + shift;  %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.2416 --> 24 actin, 16 myosin, Barezone of 24.37 nm between last myosin node and M-line
    if a_nodes == 24 && m_nodes == 16
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 48.37e-9 + shift; %59.04e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %24.37e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.3220 --> 32 actin, 20 myosin
    if a_nodes == 32 && m_nodes == 20
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 60.65e-9 + shift; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.4025 --> 1/2 HALF SARCOMERE 40 actin, 25 myosin, Barezone of 38.5 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 40 && m_nodes == 25
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 75.8e-9 + shift; %100e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %38.5e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
        
    %TEST CASE 8.4830 --> 48 actin, 30 myosin
    if a_nodes == 48 && m_nodes == 30
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 90.98e-9 + shift; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.6440 --> 1/2 HALF SARCOMERE 60 actin, 38 myosin, Barezone of 121.3 nm between last myosin node and M-line, 100 nm between last actin and M-line 
    if a_nodes == 64 && m_nodes == 40
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 121.3e-9 + shift; %100e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %38.5e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end
    
    %TEST CASE 8.8050 --> FULL HALF SARCOMERE 80 actin, 50 myosin, Barezone of 77 nm between last myosin node and M-line, 200 nm between last actin and M-line
    if a_nodes == 80 && m_nodes == 50
        N_locations(1, tot_nodes) = N_locations(1, a_nodes+1) + 151.6e-9 + shift; %137.3e-9; <- what was being used prior to 2/4/2022 %200e-9; %M-line location
        N_locations(1, (tot_nodes-1)) = N_locations(1, tot_nodes) - m0; %61.6e-9; %myosin closest to M-line location
        N_locations(1, ((a_nodes+2):(tot_nodes-2)))= N_locations(1, tot_nodes-1) -  [m0*(m_nodes-1):-m0:m0]; 
        xf = N_locations(1, tot_nodes)+xf0;
        titin0 = N_locations(1, a_nodes+2) - N_locations(1, 1); %titin connects Z-line to m1 node
        xhs = N_locations(1, tot_nodes) - N_locations(1, 1); %when titin connects Z-line to m1 node 
    end

    



%}

%% Save Constants
% filename = ['initial_geometry_',...
%             'testcase=', num2str(testcase), ',', ...
%             'xf=', num2str(m_nodes), ',', ...
%             'titin0=', num2str(a0), ',', ...
%             'xhs=', num2str(m0), '.mat'];
filename = ['ainitial_geometry_','testcase=', num2str(testcase), '_', ...
            'a_nodes=', num2str(a_nodes), ',', ...
            'm_nodes=', num2str(m_nodes), ',', ...
            'a0=', num2str(a0), ',', ...
            'm0=', num2str(m0), ',', ...
            'xb0=', num2str(xb0), ',', ...
            'xf0=', num2str(xf0), ',',  ...
            't_end=', num2str(t_end), ',', ...
            'dt=', num2str(dt), '.mat']; %shifted-0.5m0 _optimalOverlap
  
save(filename);






