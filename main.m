%% THIS IS A MINIMUM WORKING EXAMPLE

%% INCLUDE ALL THE NECESSARY FILES
    [path,~,~] = fileparts(matlab.desktop.editor.getActiveFilename) ;
    addpath(genpath(path))

%% MESH GENERATION SCRIPT
    run buildDistMesh

%% GLOBAL DIC
    run globalDIC