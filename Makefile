
unfolding_mom_v7_density_orig:
	g++ -g -pthread -m64 -Wno-deprecated -std=c++17 -I${ROOTSYS}/include -o unfolding_mom_v7_density_orig unfolding_mom_v7_density_orig.cc `root-config --cflags --libs`

clean:
	rm unfolding_mom_v7_density_orig
