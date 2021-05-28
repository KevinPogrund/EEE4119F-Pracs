%% This function will animate the robot dynamics
% 2 April 2019
% update on 17 March 2020: made compatible with Octave
%Edited by Kai Brown and Kevin Pogrund
%        BRWKAI001 and PGRKEV001
%%
close all;
clc;
clear all;
robotDynamics;
clc;
sim('prac4_Q2.slx');
disp('Kai Brown and Kevin Pogrund');
disp('BRWKAI001 and PGRKEV001');

%% Positions

% if not simming above then
% psi = out.yout.getElement('psi1');
% phi = out.yout.getElement('phi2');

%both are sampled at 10ms
psi = ans.yout.getElement('psi1');
sim_psi1 = psi.Values.Data;
phi = ans.yout.getElement('phi2');
sim_phi2 = phi.Values.Data;

%L1 = 0.5 blue
%L2 = 0.2 red
pos1 = [zeros(length(sim_psi1),1), zeros(length(sim_psi1), 1), L1_*ones(length(sim_psi1),1)];

pos2 = [-L2_*cos(sim_phi2).*sin(sim_psi1), L2_*cos(sim_phi2).*cos(sim_psi1), L2_*sin(sim_phi2)+L1_];


%% Create Figure Handles
figure
view(3);
axis([-0.8 0.8 -0.8 0.8 0 0.8])
grid on
hold on;

% The rod
h1 = line('Color', 'b', 'LineWidth', 1);
h2 = line('Color', 'r', 'LineWidth', 1);

% general plot setup
xlabel({'X Position (m)'},'FontSize',14,'FontName','AvantGarde');
ylabel({'Y Position (m)'},'FontSize',14,'FontName','AvantGarde');
zlabel({'Z Position (m)'},'FontSize',14,'FontName','AvantGarde');

title({'Robot Dynamics'},'FontWeight','bold','FontSize',20,...
    'FontName','AvantGarde');

%% Update the rod and masses in the plot in a loop
for i = 1:length(sim_phi2)
    
    set(h1, 'XData', [0, pos1(i,1)]);
    set(h1, 'YData', [0, pos1(i,2)]);
    set(h1, 'ZData', [0, pos1(i,3)]);
    
    set(h2, 'XData', [pos1(i,1), pos2(i,1)]);
    set(h2, 'YData', [pos1(i,2), pos2(i,2)]);
    set(h2, 'ZData', [pos1(i,3), pos2(i,3)]);
    
    pause(0.1);  % the time per loop is (calculation/render time) + (pause)
                 % this doesn't need to be done properly -- we'll work on
                 % that in part 2 of this prac
    drawnow 
   
    xlim('manual')
    ylim('manual')
    zlim('manual')
end
disp('done')