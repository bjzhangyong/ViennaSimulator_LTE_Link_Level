clear all
close all
warning off
clc
calling_dir = pwd;
disp('Please type in the full path to your Vienna System Level Simulator directory');
SL_dir = input('Directory path: ','s');
% SL_dir = 'D:\LTE System Level Simulator\trunk';
cd(SL_dir);
cd(calling_dir);
disp('Please type now in the full path to your Vienna Link Level Simulator directory');
LL_dir = input('Directory path: ','s');
% LL_dir = 'D:\LTE Link Level Simulator\LTE_sim_SVN';
cd(LL_dir);
cd(calling_dir);

%%
disp('Starting the Link Level Simulation');
pause(2.5)
cd([LL_dir, '\paper scripts']);
Reproducibility_Schedulers_LL;

%%
disp('Starting the System Level Simulation');
pause(2.5)
cd([SL_dir, '\paper scripts']);
Reproducibility_Schedulers_SL;

cd(calling_dir);
Reproducibility_Schedulers_plot;