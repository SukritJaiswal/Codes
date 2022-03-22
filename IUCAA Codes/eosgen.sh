gcc calccoup_isovec.c -o calccoupisovec -Wall -lgsl -lgslcblas -lm
./calccoupisovec calccoupisovec_test.out calccoup_isovec.out
gcc calceos_anm_isovec.c -o calceosanmisovec -Wall -lgsl -lgslcblas -lm
./calceosanmisovec calceosanmisovec_test.out calceosanmisovec.out calcfracanmisovec.out
