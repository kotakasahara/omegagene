========================
In/Out Files
========================

:Author: Kota Kasahara

------------------------------------
入力ファイル
------------------------------------

1. システム設定 .cfg
2. シミュレーション設定 .cfg
3. 構造情報 .pdb
4. 初期構造/速度設定 .restart
5. トポロジー設定 .tpl
6. SHAKE設定 .shk
7. V-McMD設定(2種類) .inp, .vert
8. システム情報 .cls

システム定義ファイル
------------------------------------

mdinput_generator.py 用の入力ファイルです。これらの情報を元にシステム情報ファイル .cls を生成します。
1行につき1エントリーで、各エントリーではキーワードと値のセットで設定を入力します。

* --fn-i-tpl          et1.tpl
  * トポロジー設定ファイル名を相対または絶対パスで指定します。
  * ファイルの詳細については後述します。
* --fn-i-initial-pdb  et1.pdb 
  * 構造情報を含む.pdbファイルへのパスを指定します。
  * ファイルの詳細については後述します。
* --fn-i-restart      et1.restart
  * 初期座標および速度を含むrestartファイルへのパスを指定します。
  * ファイルの詳細については後述します。
* --cell-x            61.2425 
* --cell-y            61.2425
* --cell-z            61.2425
  * 周期境界セルのX,Y,Z軸のセル長をangestrome単位で指定します。
* --fn-i-shake        system.s
  * SHAKE設定ファイルへのパスを指定します。
* --fn-i-ttp-v-mcmd-inp      ttp_v_mcmd.inp
* --fn-i-ttp-v-mcmd-initial  start.vert    
  * V-McMDファイルへのパスを指定します。

シミュレーション設定ファイル
------------------------------------

Celeste用入力ファイルです。このファイルと、システム情報ファイル .cls をセットでCelesteに入力します。

md_input.cfg::

* --mode                  md
  * 現状 "md" のみを許可

* --integrator            leapfrog-presto
  * 現状下記2つのみ許可
    * leapfrog-presto
      * Psygeneと同様の積分アルゴリズム
      * 熱浴なし、"scaling"熱浴、SHAKEに対応
    * zhang
      * Zhang [Zhang1997]_ の積分アルゴリズム
      * "hoover-evans"熱浴のみ可。SHAKE不可。
* --thermostat            scaling
  * 現状下記3つのみ許可
    * none
      * 熱浴無し。マイクロカノニカルアンサンブル
    * scaling
      * 温度スケーリング
    * hoover-evans
      * Hoover-Evans熱浴。--integratorはzhangのみ対応。
* --cutoff                12.0
  * 非結合性ペアポテンシャルのカットオフ距離
* --n-steps               10 
  * 計算ステップ数
* --time-step             2.0
  * 積分時間
* --electrostatic         zero-dipole
  * 現状 "zero-dipole"のみ可
* --ele-alpha             0.0
  * Zero-multipole summation法のダンピング係数
  * 現状、GPU版は0.0のみ対応
* --temperature           300        
  * 熱浴使用時の温度
* --print-interval-log    1          
  * ログの出力間隔
* --print-interval-coord  1          
  * 座標情報の出力間隔
* --fn-o-coord            et1.trr    
  * 座標情報の出力先ファイル名
* --format-o-coord        gromacs    
  * 座標情報の出力形式
  * gromacs
  * presto
* --fn-o-log              et1.log    
  * 未実装
* --fn-o-energy           et1.ene    
  * 未実装
* --nsgrid-cutoff         13.0       
  * 近傍探索用カットオフ
* --nsgrid-update-intvl   50
  * 近傍探索の間隔
* --expanded-ensemble             v-mcmd
  * none
    * 拡張アンサンブルを使用しない
  * v-mcmd
    * TTP-V-McMD法[Higo2013]_ 
* --fn-o-vmcmd-log                ttp_v_mcmd.out
  * V-McMDの仮想状態トラジェクトリの出力先ファイル名
* --fn-o-expand-lambda            mule.ene       
  * 拡張アンサンブルの反応座標トラジェクトリの出力先ファイル名
* --print-interval-expand-lambda  1              
  * 反応座標トラジェクトリの出力間隔
* --format-o-expand-lambda        ascii          
  * 反応座標トラジェクトリの出力形式
  * ascii
  * binary

.. [Zhang1997]_ Zhang "Operator-splitting integrators for constant-temperature molecular dynamics" J. Chem. Phys. 106 (14) 1997


初期構造/速度設定ファイル
------------------------------------
myPresto/PsygeneのRestart fileに準拠。
附属の presto_generate_velocities.py スクリプトでも生成出来る。

トポロジー設定ファイル
------------------------------------
myPresto/Psygeneに準拠。
myPrestoパッケージ附属のtplgeneプログラムで生成出来る。

SHAKE設定ファイル
------------------------------------
myPresto/Psygeneに準拠。
myPrestoパッケージ附属のSHAKEinpプログラムで生成出来る。

V-McMD設定ファイル
------------------------------------
myPresto/Psygeneに準拠。
* --fn-i-ttp-v-mcmd-inp      ttp_v_mcmd.inp
* --fn-i-ttp-v-mcmd-initial  start.vert    

------------------------------------
出力ファイル
------------------------------------

1. 構造トラジェクトリファイル
2. 最終構造/速度ファイル
3. V-McMD出力ファイル(2種類)

