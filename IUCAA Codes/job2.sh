gcc calccoup_isovec_1.c -o calccoupisovec -Wall -lgsl -lgslcblas -lm
./calccoupisovec calccoupisovec_test.out calccoup_isovec.out
gcc calceos_anm_isovec_1.c -o calceosanmisovec -Wall -lgsl -lgslcblas -lm
./calceosanmisovec calceosanmisovec_test.out calceosanmisovec.out calcfracanmisovec.out
#python3 crust_core_matching.py
#gcc tov_phi.c -o tovphi -Wall  -lm
#./tovphi tovphi_test.out  tovphi_plot.out
#gcc pmode.c -o pmode -Wall -lm
#./pmode pmode.out
