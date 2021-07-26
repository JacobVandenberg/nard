function [filename] = write_config(conf,filename, force)
%WRITE_CONFIG writes a nard config to an h5 file
%   Detailed explanation goes here

if nargin < 3
    force = false;
end

if isfile(filename)
    if force
        
    else
        if input(sprintf("Config file (%s) already exists! Overwrite? (y/n): ", filename), 's') ~= 'y'
            fprintf("Aborting.")
            error("File Already Exists.")
        end
    end
    delete(filename)
end


h5create(filename, "/x", numel(conf.x), 'Datatype', 'double');
h5write(filename, "/x", conf.x)

h5create(filename, "/y", numel(conf.y), 'Datatype', 'double');
h5write(filename, "/y", conf.y)

h5create(filename, "/diffusion_consts", numel(conf.diffusion_consts), 'Datatype', 'double');
h5write(filename, "/diffusion_consts", conf.diffusion_consts)

h5create(filename, "/DBCx_plus", numel(conf.DBCx_plus), 'Datatype', 'double');
h5write(filename, "/DBCx_plus", conf.DBCx_plus)

h5create(filename, "/DBCx_minus", numel(conf.DBCx_minus), 'Datatype', 'double');
h5write(filename, "/DBCx_minus", conf.DBCx_minus)

h5create(filename, "/DBCy_plus", numel(conf.DBCy_plus), 'Datatype', 'double');
h5write(filename, "/DBCy_plus", conf.DBCy_plus)

h5create(filename, "/DBCy_minus", numel(conf.DBCy_minus), 'Datatype', 'double');
h5write(filename, "/DBCy_minus", conf.DBCy_minus)

h5create(filename, "/rparams", numel(conf.rparams), 'Datatype', 'double');
h5write(filename, "/rparams", conf.rparams)

h5create(filename, "/user_params", numel(conf.user_params), 'Datatype', 'double');
h5write(filename, "/user_params", conf.user_params)

h5create(filename, "/IC", size(conf.IC), 'Datatype', 'double');
h5write(filename, "/IC", conf.IC)

h5create(filename, "/iparams", numel(conf.iparams), 'Datatype', 'int64');
h5write(filename, "/iparams", conf.iparams)

h5create(filename, "/DBCx_plus_mask", numel(conf.DBCx_plus_mask), 'Datatype', 'int64');
h5write(filename, "/DBCx_plus_mask", conf.DBCx_plus_mask)

h5create(filename, "/DBCx_minus_mask", numel(conf.DBCx_minus_mask), 'Datatype', 'int64');
h5write(filename, "/DBCx_minus_mask", conf.DBCx_minus_mask)

h5create(filename, "/DBCy_plus_mask", numel(conf.DBCy_plus_mask), 'Datatype', 'int64');
h5write(filename, "/DBCy_plus_mask", conf.DBCy_plus_mask)

h5create(filename, "/DBCy_minus_mask", numel(conf.DBCy_minus_mask), 'Datatype', 'int64');
h5write(filename, "/DBCy_minus_mask", conf.DBCy_minus_mask)

h5create(filename, "/savefilename", 1, 'Datatype', 'string');
h5write(filename, "/savefilename", conf.savefilename)

h5create(filename, "/plotfilename", 1, 'Datatype', 'string');
h5write(filename, "/plotfilename", conf.plotfilename)

end

