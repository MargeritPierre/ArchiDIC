%% THIS IS A MINIMUM WORKING EXAMPLE
    clear all
    close all
    clc

%% INCLUDE ALL THE NECESSARY FILES
    [path,~,~] = fileparts(matlab.desktop.editor.getActiveFilename) ;
    addpath(genpath(path))

%% MESH GENERATION SCRIPT
    run buildDistMesh

%% MICROSCOPIC MEASUREMENTS
    run globalDIC

%% MACROSCOPIC MEASUREMENTS
    run localDIC