=========================
|Celeste| Project Roadmap
=========================

.. contents::

This page describes technical ideas for improvements of the calculation speed.
The ideas are suggestions for Goto-san of TITECH.


------------------------------
OpenMP Implementation on SHAKE
------------------------------

[SHAKEをOpenMP化する]

研究としてあまり楽しくはないかも知れませんが、簡単に出来てかつ効果がありそうなところです。

SHAKEは、

``ConstraintShake.cpp : ConstraintShake::apply_constraint()``

関数を見て頂くと分かりますが、
SHAKE単位毎（CH3とか、CH2とか、重原子と水素原子のセット）に単純にループを回していっこずつ処理しています。
ここは割と簡単にOpenMPに載せられると思いますので、楽してポイントを稼げる所だと思います。


------------------------------------------------------------------
Parallelization Neighbor Search and Potential Calculation Routines
------------------------------------------------------------------

[Neighbor Search とポテンシャル計算を並列化させる]

今のところ、数ステップおきにセルに切ってペアを作る処理を、ポテンシャル計算の前にやっています。
しかしポテンシャル計算をGPUでやりながら、1step前の座標を使ってセルペア列挙をCPUで行えるはずです。
この場合、F_ECPモードのルーチンを使うことになります。
またここでは、セル1に近いセルを列挙、セル2に近いセルを列挙…と逐次的にループで処理しているので、ここも比較的簡単にOpenMPに載せられると思います。


--------------------------------------------------------------
Modifying data structures for effective data access on SIMD(T)
--------------------------------------------------------------

[データ構造]

いまのところ、一本の配列に

>  [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]

とNこの原子のXYZ座標を押し込んでいます。
またこれとは別に各原子の電荷パラメータを要素数Nの一本の配列に入れてます。

>  [q1, q2, ..., qN]

この値はシミュレーションを通して一定です。
これをGPU上でfloat4 の配列に押し込んで、

>  [ float4(x1, y1, z1, q1), float4(x2, y2, z2, q2), ..., float4(xN,yN,zN,qN) ]

としてから計算をしているのですが、いまいちうまくない気がしています。

たとえば、

>  [x1, x2, x3, ..., xN]
>  [y1, y2, y3, ..., yN]
>  [z1, z2, z3, ..., zN]
>  [q1, q2, q3, ..., zN]

という四本の配列にしておいて、
セル内の8原子分のデータを一気に読んでshared memoryに入れておけばわざわざfloat4に入れなくても良いような気がします。

またここはCPUのSIMD化にも必要ですので、そのうちCPU部分のデータ構造も書き換えたいと思っています。

ver.0.37.a.1 (v037a_arraymod ブランチ）にCPU部分のデータ構造を修正したバージョンを置いておきました。


--------------------------------------------------------
Function Calls for ZM in Pair-wise Potential Calculation
--------------------------------------------------------

[二体ポテンシャル計算におけるZMの取り扱い]

