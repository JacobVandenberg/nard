function [] = nard2(reaction_name, config_file)
%NARD2 Summary of this function goes here
%   Detailed explanation goes here
    randomstr = char(randi([65 65+25],1,20));
    filename = sprintf(".nard2_tempconfig_%s.h5", randomstr);
    written_file = false;
    if (class(config_file) == "struct")
        written_file = true;
        write_config(config_file, filename)
        config_file = filename;
    end
    try
        system(sprintf("nard2 run %s %s", reaction_name, config_file))
    end
    if written_file
        delete(filename);
    end
end

