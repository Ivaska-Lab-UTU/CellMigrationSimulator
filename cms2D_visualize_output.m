%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aleksi Isomursu
% University of Turku
% February 2023

% Script for visualizing the output of the Odde lab 2D Cell Migration
% Simulator (Bangasser et al., 2017), produced by cms2Dv1_0.m and
% analyzeCMSoutputv1_0.m.

% Import the MATLAB data representing a single simulated cell (e.g.
% cms2Dv1_0_Ksub0.5_Nmotors7500_Nclutch7500_run1.m), set the desired
% temporal resolution (relative to the original data, e.g. original tint 
% of 1 sec and interval of 60 means that one data point per a minute of
% simulation will be visualized) and size of the resulting image(s), 
% and run. Note: all coordinates are converted from nm to um.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear; close all; clc
%%
% PARAMETERS and INPUT
% Define the temporal sampling resolution, in terms of the original 'tint'.
interval = 60;

% Define image dimensions (for x and y, in um).
imgsize = 300;

% Define the nuclear radius (in um)
nucrad = 10;

% Define the .mat file to visualize.
[inputfile,inputpath] = uigetfile('*.mat', ['Select a CMS .mat file to ' ...
    'visualize.']);
[path,filename,ext] = fileparts(strcat(inputpath,inputfile));
% Define (and create) the output folder.
outputdir = strcat(path,'\Visualizations\',filename,'\');
if ~exist(outputdir, 'dir')
    mkdir(outputdir);
end

%%
% Load the data.
load(strcat(inputpath, inputfile));

imgindex = 1;
% Loop through the selected simulation time points to analyze them.
for i = 1:interval:length(cellpos)
    % Record the position of the cell body at the given time point.
    coordcell = cellpos(:,i)/1000;
    % Initialize an empty array for the module coordinates at the given
    % time point.
    coordmodules = zeros(2,modnum(i));
    % Record the positions of the modules, taking into account the
    % substrate displacement.
    for j = 1:modnum(i)
        coordmodules(:,j) = module(j).xsub(:,i)/1000;
    end

    % Start drawing the figure.
    f = figure('Color', [0.95 0.95 0.95]);
    f.WindowState = 'maximized';            % full screen
    daspect([1 1 1]);                       % image aspect ratio
    xlim([-imgsize/2, imgsize/2]);          % axes
    ylim([-imgsize/2, imgsize/2]);
    hold on;

    % Draw the nucleus.
    theta = 0:pi/30:2*pi;
    xnuc = nucrad * cos(theta) + coordcell(1,1);        % r,x
    ynuc = nucrad * sin(theta) + coordcell(2,1);        % r,y
    fill(xnuc, ynuc, [0.3 0.3 0.3],'LineStyle','none');

    % Add the individual modules.
    for j = 1:length(coordmodules)
        % Resolve the points where the protrusion outlines contact the
        % sides of the nucleus (tangents).
        p1 = [coordcell(1,1), coordcell(2,1)];          % nucleus
        p2 = [coordmodules(1,j), coordmodules(2,j)];    % module

        % a,b~x,y representing the points of intersection
        syms a b
        q = [a,b];
        % 'q' lies on the circle (equation of the circle).
        eqn1 = (q(1)-p1(1))^2 + (q(2)-p1(2))^2 - nucrad^2 == 0;
        % Pythagoras applies (tangent line and the perpendicular slope ~r).
        eqn2 = (p2(2)-q(2))/(p2(1)-q(1)) == -(q(1)-p1(1))/(q(2)-p1(2));

        res = solve(eqn1,eqn2);
        xmod = double(res.a);
        ymod = double(res.b);
        % If q could not be resolved (i.e. module under the nucleus), skip.
        if ~isreal(xmod) || ~isreal(ymod)
            continue;
        end

        % Draw the module based on the coordinates of its three vertices.
        % Arguments: ([x1,x2,x3],[y1,y2,y3],color,...)
        fill([p2(1) xmod(1) xmod(2)],[p2(2) ymod(1) ymod(2)], ...
            [0.3 0.3 0.3],'LineStyle','none');
    end

    % Save the results in the output folder and pad the index with zeros.
    fr = getframe;
    imwrite(fr.cdata,strcat(outputdir,'img',num2str(imgindex,'%04.f'), ...
        '.tiff'));
    hold off;
    close(f);
    imgindex = imgindex + 1;
end