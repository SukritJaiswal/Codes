gcc calccoup_isovec.c -o calccoupisovec -Wall -lgsl -lgslcblas -lm
./calccoupisovec calccoupisovec_test.out calccoup_isovec.out
gcc calceos_anm_isovec.c -o calceosanmisovec -Wall -lgsl -lgslcblas -lm
./calceosanmisovec calceosanmisovec_test.out calceosanmisovec.out calcfracanmisovec.out
python3 crust_core_matching.py
gcc tov_phi.c -o tovphi -Wall  -lm
./tovphi tovphi_test.out  tovphi_plot.out
#gcc modes_corr.c -o loopsol -Wall -lgsl -lgslcblas -lm
#./loopsol loopsol_test.out  loopsol_plot.out
gcc modes_lowm.c -o lowm -Wall -lm
./lowm pmode_lowmass.out
#gcc pmode.c -o pmode -Wall -lm
#./pmode pmode.out
