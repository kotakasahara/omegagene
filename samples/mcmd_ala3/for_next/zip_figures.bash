#!/bin/bash

head -1 current_situation > aaa
mdn=`more aaa`
rm aaa
head -2 current_situation > aaa
tail -1 aaa > bbb
nvs=`more bbb`
rm aaa bbb
if [ $mdn -lt 10 ]; then
  mdn_str="0${mdn}"
else
 mdn_str="${mdn}"
fi

tar cvzf fig.tar.gz  v_distrib/cal${mdn_str}_mcmd1/v_distrib.png  fit_pmc_entire/fitpmc*.png
