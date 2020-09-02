#

  rm toto/*
  @ ii = 1
  while( $ii <= 12)
    echo $ii
    cp ${ii}/dden.dat  toto/m${ii}
    @ ii ++
  end

  tar cvf toto.tar toto
  cp toto.tar ../../../
 
