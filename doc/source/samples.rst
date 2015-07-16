========================
Samples
========================

:Author: Kota Kasahara

--------------------------
サンプルファイル
--------------------------

サンプル計算ファイルとしてEndothelin-1 dimer (PDB-ID: 1t7h) を付属しています。

:: 

  python2.7 ${CELESTETK}/mdinput_generator.py -i gen_input.cfg -o et1.cls
  celeste_gpu --inp et1.cls --cfg md_input.cfg

