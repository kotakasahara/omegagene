set style data line
set key outside
set size 1,1
set origin 0,0
plot [4.7:20.0] []    's1/e1.pdf' title 's1', \
                      's2/e1.pdf' title 's2', \
                      's3/e1.pdf' title 's3', \
                      's4/e1.pdf' title 's4', \
                      's5/e1.pdf' title 's5', \
                      's6/e1.pdf' title 's6', \
                      's7/e1.pdf' title 's7'
