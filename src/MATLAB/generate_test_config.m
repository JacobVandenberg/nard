conf.x = linspace(0, 1, 5);
conf.y = conf.x + 1;
conf.diffusion_consts = [1.0 30.0];

conf.DBCx_plus = [0.0 0.0];
conf.DBCx_minus = [0.0 0.0];
conf.DBCy_plus = [0.0 0.0];
conf.DBCy_minus = [0.0 0.0];

conf.rparams = zeros(64, 1);

conf.rparams(1:3) = [0.00001 1.0 30.0];
conf.iparams = zeros(64, 1);
conf.iparams(1:5) = [1000 100000000 0 0 2];
conf.iparams = int64(conf.iparams);

conf.DBCx_plus_mask = int64([0 0]);
conf.DBCx_minus_mask = int64([0 0]);
conf.DBCy_plus_mask = int64([0 0]);
conf.DBCy_minus_mask = int64([0 0]);

conf.savefilename = "src/tests/test_config_result.h5";
conf.plotfilename = "src/tests/test_config_plot.png";

h5create("src/tests/test_config.h5", "/x", numel(conf.x), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/x", conf.x)

h5create("src/tests/test_config.h5", "/y", numel(conf.y), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/y", conf.y)

h5create("src/tests/test_config.h5", "/diffusion_consts", numel(conf.diffusion_consts), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/diffusion_consts", conf.diffusion_consts)

h5create("src/tests/test_config.h5", "/DBCx_plus", numel(conf.DBCx_plus), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/DBCx_plus", conf.DBCx_plus)

h5create("src/tests/test_config.h5", "/DBCx_minus", numel(conf.DBCx_minus), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/DBCx_minus", conf.DBCx_minus)

h5create("src/tests/test_config.h5", "/DBCy_plus", numel(conf.DBCy_plus), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/DBCy_plus", conf.DBCy_plus)

h5create("src/tests/test_config.h5", "/DBCy_minus", numel(conf.DBCy_minus), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/DBCy_minus", conf.DBCy_minus)

h5create("src/tests/test_config.h5", "/rparams", numel(conf.rparams), 'Datatype', 'double');
h5write("src/tests/test_config.h5", "/rparams", conf.rparams)

h5create("src/tests/test_config.h5", "/iparams", numel(conf.iparams), 'Datatype', 'int64');
h5write("src/tests/test_config.h5", "/iparams", conf.iparams)

h5create("src/tests/test_config.h5", "/DBCx_plus_mask", numel(conf.DBCx_plus_mask), 'Datatype', 'int64');
h5write("src/tests/test_config.h5", "/DBCx_plus_mask", conf.DBCx_plus_mask)

h5create("src/tests/test_config.h5", "/DBCx_minus_mask", numel(conf.DBCx_minus_mask), 'Datatype', 'int64');
h5write("src/tests/test_config.h5", "/DBCx_minus_mask", conf.DBCx_minus_mask)

h5create("src/tests/test_config.h5", "/DBCy_plus_mask", numel(conf.DBCy_plus_mask), 'Datatype', 'int64');
h5write("src/tests/test_config.h5", "/DBCy_plus_mask", conf.DBCy_plus_mask)

h5create("src/tests/test_config.h5", "/DBCy_minus_mask", numel(conf.DBCy_minus_mask), 'Datatype', 'int64');
h5write("src/tests/test_config.h5", "/DBCy_minus_mask", conf.DBCy_minus_mask)

h5create("src/tests/test_config.h5", "/savefilename", 1, 'Datatype', 'string');
h5write("src/tests/test_config.h5", "/savefilename", conf.savefilename)

h5create("src/tests/test_config.h5", "/plotfilename", 1, 'Datatype', 'string');
h5write("src/tests/test_config.h5", "/plotfilename", conf.plotfilename)

