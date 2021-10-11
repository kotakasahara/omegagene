
#g++ prob_dblwell_3d.cpp -o prob_dblwell_3d
#./prob_dblwell_3d 10 > prob_dblwell_3d_10.dat
./prob_dblwell_3d 11 > prob_dblwell_3d_11.dat
./prob_dblwell_3d 12 > prob_dblwell_3d_12.dat
./prob_dblwell_3d 13 > prob_dblwell_3d_13.dat
./prob_dblwell_3d 14 > prob_dblwell_3d_14.dat
#./prob_dblwell_3d 15 > prob_dblwell_3d_15.dat
#./prob_dblwell_3d 20 > prob_dblwell_3d_20.dat
#./prob_dblwell_3d 25 > prob_dblwell_3d_25.dat
#./prob_dblwell_3d 30 > prob_dblwell_3d_30.dat


#mv prob_dblwell_3d_10.dat tmp.dat
#cat vcmd.inp tmp.dat > prob_dblwell_3d_10.dat
mv prob_dblwell_3d_11.dat tmp.dat
cat vcmd.inp tmp.dat > prob_dblwell_3d_11.dat
mv prob_dblwell_3d_12.dat tmp.dat
cat vcmd.inp tmp.dat > prob_dblwell_3d_12.dat
mv prob_dblwell_3d_13.dat tmp.dat
cat vcmd.inp tmp.dat > prob_dblwell_3d_13.dat
mv prob_dblwell_3d_14.dat tmp.dat
cat vcmd.inp tmp.dat > prob_dblwell_3d_14.dat
#mv prob_dblwell_3d_15.dat tmp.dat
#cat vcmd.inp tmp.dat > prob_dblwell_3d_15.dat
#mv prob_dblwell_3d_20.dat tmp.dat
#cat vcmd.inp tmp.dat > prob_dblwell_3d_20.dat
#mv prob_dblwell_3d_25.dat tmp.dat
#cat vcmd.inp tmp.dat > prob_dblwell_3d_25.dat
#mv prob_dblwell_3d_30.dat tmp.dat
#cat vcmd.inp tmp.dat > prob_dblwell_3d_30.dat
