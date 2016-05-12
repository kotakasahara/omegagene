====================================
二体力計算ルーチンに関するメモ
====================================

:Author: Kota Kasahara

------------------------------------
二体力計算に至るまでのながれ
------------------------------------

Celesteでは原子の座標をまずMmSystemクラスのcrd配列に読み込みます。
MmSystem.crd配列は [原子数, 3] の二次元配列で、各原子のx,y,z座標を記録します。::

  crd[atomid][0] = x;
  crd[atomid][1] = y;
  crd[atomid][2] = z;

計算の前に、crdの内容はSubBoxクラスのcrdへと渡されます。
これは将来的にはプロセス並列を実装するためのもので、分割した空間の一部のみがSubBoxへ渡されることになります。
現バージョンでは空間分割を実装していないので、MmSystem.crdとSubBox.crdは同じ内容になります。
但しSubBox.crdは長さが(原子数*3)の一次元配列になっており、xyzxyzxyz...の順番で値が記録されます。
またSubBox.atomids配列にはMmSystem側との原子IDの対応を記録します。::

  SubBox.atomids[SubBoxでの原子ID] = MmSystemでの原子ID
  MmSystem.crd[SubBox.atomids[i]][0] = crd[i*3+0]

------------------------------------------
Neighbor searchなしでの計算ルーチンの概要
------------------------------------------

Neighbor search を行わずに二体力計算を行う場合には、crd配列に対して二重ループを回すことで全対全の二体力を計算します。

また特定の原子ペアについては二体力計算をスキップしなくてはならないため、条件分岐を行います。（後述）

具体的な二体力の計算式はForceFieldクラスのcalc_pairwise関数に座標とパラメータを渡すことで実行します。
原子間距離がカットオフ以内だった場合は戻り値0が得られますので、そのときに限り計算結果を記録します。

SubBox\:\:calc_energy_pairwise_wo_neighborsearch::

  int SubBox::calc_energy_pairwise_wo_neighborsearch(){
    for(int atomid1 = 0, atomid1_3=0; atomid1 < n_atoms_exbox; atomid1++, atomid1_3+=3){
      real crd1[3] = {crd[atomid1_3], crd[atomid1_3+1], crd[atomid1_3+2]};
      for(int atomid2 = 0, atomid2_3=0; atomid2 < atomid1; atomid2++, atomid2_3+=3){
        real crd2[3] = {crd[atomid2_3], crd[atomid2_3+1], crd[atomid2_3+2]};
        bool flg=true;
        for(int i = atomid1 * max_n_nb15off;
            i < atomid1 * max_n_nb15off + max_n_nb15off;
            i++){
          if(nb15off[i] == atomid2) flg = false; 
        }
        if(!flg) continue;

        real_pw tmp_ene_vdw = 0.0;
        real_pw tmp_ene_ele = 0.0;
        real_fc tmp_work[3] = {0.0, 0.0, 0.0};
        real param_6term  = lj_6term[atom_type[atomid1]  * n_lj_types + atom_type[atomid2]];
        real param_12term = lj_12term[atom_type[atomid1] * n_lj_types + atom_type[atomid2]];

        for(int d=0; d < 3; d++){
          if(crd2[d]-crd1[d] >= pbc->L_half[d])
            crd2[d] -= pbc->L[d];
          else if(crd2[d]-crd1[d] <= -pbc->L_half[d])
            crd2[d] += pbc->L[d];
        }
        if(ff.calc_pairwise(tmp_ene_vdw, tmp_ene_ele, tmp_work,
                            crd1, crd2,
                            param_6term, param_12term,
                            charge[atomid1],
                            charge[atomid2])==0){
          pote_vdw += tmp_ene_vdw;
          pote_ele += tmp_ene_ele;
          work[atomid1_3] += tmp_work[0];
          work[atomid1_3+1] += tmp_work[1];
          work[atomid1_3+2] += tmp_work[2];
          work[atomid2_3] -= tmp_work[0];
          work[atomid2_3+1] -= tmp_work[1];
          work[atomid2_3+2] -= tmp_work[2];
        }
      }
    }
    return 0;
  }


------------------------------------
結合原子ペアの除外
------------------------------------

結合性ポテンシャルが定義されている原子ペアについては、非結合性の二体ポテンシャルを計算しません。
つまり直接共有結合している隣の原子や、共有結合を2つおよび3つ挟んだ近傍の原子ペアはスキップします。
nb15off配列にどの原子ペアでの計算をスキップするかを記録してあるので、これを随時参照して判定します。

nb15offは長さが (原子数 * max_n_nb15off) のint型配列です。
配列のi*max_n_nb15off 番目からmax_n_nb15off個の値は、i番目の原子との相互作用を計算しない原子IDのリストになっています。
i-jペアに対して判定を行う際は、 この範囲にjが含まれていればスキップを行います。

------------------------------------
二体力の計算
------------------------------------

計算には原子座標の他にパラメータcharge、atom_type、lj_6term、lj_12termを利用します。
chargeは各原子の部分電荷を示すfloat配列です。
atom_typeは各原子がどのタイプに属するかを指定します。原子はn_lj_types種類（40程度）の原子タイプに分かれ、原子タイプのペアごとにvan der Waalsパラメータが決まります。lj_6termとlj_12termはペア毎のパラメータを示す配列で、それぞれ (n_lj_types * n_lj_types) の長さを持っています。

::

  real param_6term  = lj_6term[atom_type[atomid1]  * n_lj_types + atom_type[atomid2]];
  real param_12term = lj_12term[atom_type[atomid1] * n_lj_types + atom_type[atomid2]];

なお静電相互作用の計算ルーチンはZeroMultipoleSum.cppに実装されています。ZeroDipoleSum.cppは現状なにもしていません。

