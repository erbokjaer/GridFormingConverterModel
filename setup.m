clear

% Add plotFunctions directory to the MATLAB path
addpath('plotFunctions');
addpath('helperFunctions')

% Initialize random seed for reproducibility
rng(203829);

% Set figure plasment
set(groot, 'DefaultFigureWindowStyle', 'docked');
% set(groot, 'DefaultFigureWindowStyle', 'normal');
% set(groot, 'DefaultFigureWindowStyle', 'modal');

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


% Create an output folder if it does not exist
export = 0;

outputFolder = 'output';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

fig_type = ".png";

%
set(groot, ...
    'DefaultAxesFontSize', 12, ...      
    'DefaultTextFontSize', 12, ...      
    'DefaultLineLineWidth', 1.5);      

set(groot, ...
    'DefaultAxesLineWidth', 1.2, ...      
    'DefaultAxesGridLineStyle','-', ...
    'DefaultAxesGridAlpha',0.3 );


function C = uni_colors()

C.corporate_red = [153 0 0]/255;
C.blue          = [47 62 234]/255;
C.bright_green  = [31 208 130]/255;
C.navy_blue     = [3 15 79]/255;
C.yellow        = [246 208 77]/255;
C.orange        = [252 118 52]/255;
C.pink          = [247 187 177]/255;
C.grey          = [218 218 218]/255;
C.red           = [232 63 72]/255;
C.green         = [0 136 53]/255;
C.purple        = [121 35 142]/255;

end

C = uni_colors;

set(groot,'defaultAxesColorOrder',[
    C.blue
    C.red
    C.green
    C.orange
    C.purple
    C.yellow
    C.navy_blue
]);

