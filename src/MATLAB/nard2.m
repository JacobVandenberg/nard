function [] = nard2(reaction_name, config_file)
%NARD2 Summary of this function goes here
%   Detailed explanation goes here
    system(sprintf("nard2 run %s %s", reaction_name, config_file))
end

