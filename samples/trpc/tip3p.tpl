TPL> TITLE
TIP3P

TPL> MOLECULES
;NUMBER OF MOLECULES :  NUM_OF_WATER
WAT                     NUM_OF_WATER    ; set number of water molecules

TPL> ATOMS
WAT
; NUMBER OF ATOMS =          3
O       OW     18 WAT         1   16.0000    1.7683   -0.8340  2  0  0 -> ;    1
      1      2                                                         ->
      0      0      0      0       0.0000    0.0000    0.0000
H1      HW     14 WAT         1    1.0080    0.0000    0.4170  1  0  0 -> ;    2
      1                                                                ->
     -1      0      0      0       0.9572    0.0000    0.0000
H2      HW     14 WAT         1    1.0080    0.0000    0.4170  0  0  0 -> ;    3
     -2     -1      0      0       0.9572  104.5200    0.0000

TPL> BONDS
WAT
; NUMBER OF BONDS =          3
         1         2    553.0000000      0.9572000 ;          1
         1         3    553.0000000      0.9572000 ;          2
         2         3    553.0000000      1.5136000 ;          3

TPL> FUNCTIONS
;NONBONDED-POTENTIALS-FOR-PREPARATION
     1     4      LENNARD-JONES-AMBER                         ;(R*-AND-E)
; TYPE OF ATOMS
;    c,c*,ca,cb,cc,cd,ck,cm,cn,cq,cr,cv,cw,cx,cy,cz,c1,c2    :   1
;    c3,4c1                                                  :   2
;    h                                                       :   3
;    ho                                                      :   4
;    hs                                                      :   5
;    hc,1h1,hn                                               :   6
;    h1                                                      :   7
;    h2                                                      :   8
;    h3                                                      :   9
;    hp                                                      :  10
;    ha,HZ                                                   :  11
;    h4                                                      :  12
;    h5                                                      :  13
;    hw                                                      :  14
;    n,n*,n2,n3,na,nB,nc,np,no,ny,1n1,2n1,2n2,3n1,3n2,4n1,nt,n1,n4,nb :  15
;    o                                                       :  16
;    o2,1o1                                                  :  17
;    ow                                                      :  18
;    oh                                                      :  19
;    os,2o1                                                  :  20
;    p,4p1,p3,p4,p5                                          :  21
;    s,2s1,ss,s2,s4,s6                                       :  22
;    sh                                                      :  23
;    im                                                      :  24
;    ip                                                      :  25
;    k                                                       :  26
;    li                                                      :  27
;    na                                                      :  28
;    rb                                                      :  29
;    cs                                                      :  30
;    mg                                                      :  31
;    c0                                                      :  32
;    zn                                                      :  33
;    f                                                       :  34
;    cl                                                      :  35
;    br                                                      :  36
;    i                                                       :  37
;    ib                                                      :  38
;    lp                                                      :  39

TPL> NONBONDS
;NUMBER OF TYPE=    39
     1     0  1   1.90800   0.086000   0.8333333   0.500; c
     2     0  1   1.90800   0.109400   0.8333333   0.500; c3
     3     0  1   0.60000   0.015700   0.8333333   0.500; h
     4     0  1   0.00000   0.000000   0.8333333   0.500; ho
     5     0  1   0.60000   0.015700   0.8333333   0.500; hs
     6     0  1   1.48700   0.015700   0.8333333   0.500; hc
     7     0  1   1.38700   0.015700   0.8333333   0.500; h1
     8     0  1   1.28700   0.015700   0.8333333   0.500; h2
     9     0  1   1.18700   0.015700   0.8333333   0.500; h3
    10     0  1   1.10000   0.015700   0.8333333   0.500; hp
    11     0  1   1.45900   0.015000   0.8333333   0.500; ha
    12     0  1   1.40900   0.015000   0.8333333   0.500; h4
    13     0  1   1.35900   0.015000   0.8333333   0.500; h5
    14     0  1   0.00000   0.000000   0.8333333   0.500; hw
    15     0  1   1.82400   0.170000   0.8333333   0.500; n
    16     0  1   1.66120   0.210000   0.8333333   0.500; o
    17     0  1   1.66120   0.210000   0.8333333   0.500; o2
    18     0  1   1.76830   0.152000   0.8333333   0.500; ow
    19     0  1   1.72100   0.210400   0.8333333   0.500; oh
    20     0  1   1.68370   0.170000   0.8333333   0.500; os
    21     0  1   2.10000   0.200000   0.8333333   0.500; p
    22     0  1   2.00000   0.250000   0.8333333   0.500; s
    23     0  1   2.00000   0.250000   0.8333333   0.500; sh
    24     0  1   2.47000   0.100000   0.8333333   0.500; IM
    25     0  1   1.86800   0.002770   0.8333333   0.500; Ip 
    26     0  1   2.65800   0.000328   0.8333333   0.500; K  
    27     0  1   1.13700   0.0183     0.8333333   0.500; Li 
    28     0  1   1.86800   0.00277    0.8333333   0.500; na 
    29     0  1   2.95600   0.00017    0.8333333   0.500; Rb 
    30     0  1   3.39500   0.0000806  0.8333333   0.500; cs 
    31     0  1   0.79260   0.8947     0.8333333   0.500; MG 
    32     0  1   1.7131    0.459789   0.8333333   0.500; ca+2  
    33     0  1   1.1000    0.012500   0.8333333   0.500; Zn+2  
    34     0  1   1.75      0.061      0.8333333   0.500; F     
    35     0  1   1.948     0.265      0.8333333   0.500; cl    
    36     0  1   2.22      0.320      0.8333333   0.500; Br    
    37     0  1   2.35      0.40       0.8333333   0.500; I     
    38     0  1   5.0       0.1        0.8333333   0.500; IB 
    39     0  1   0.0       0.000      0.8333333   0.500; Lp    

