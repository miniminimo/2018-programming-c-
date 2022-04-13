# 2018-programming-c-
c++ files
Dosing-投藥分析 
如果未投藥:每個網格 k 在第 t 期內，會發生病毒的「成長」、「 擴散」、「原生」等現象， 
(1) 「成長」： 我們假設期初時該網格已承接其上一期期末的病毒量 Q0[i]，而這些量經過一期之後
會有 s 比例的「成長量」，亦即這些承接自上期的本格總病毒量為 Q1[i]=Q0[i]*(1+s*0.01) 
(2) 「擴散」： 同理 Q0[i] 經過一期之後亦會有 q 比例擴散至與其曼哈頓距離=1 的相鄰網格 k’，因
此經過一期之後將承接自鄰近格擴散進來的總病毒量為 Q2[i]=加總的 Q0[i]*q 
(3) 「原生」： 在本期間也可能依已知的機率「原生」出 Q3 單位的原生病毒量 
綜上所述，在未投藥的情況下，網格 k 第 t 期期末總病毒量 Q0[i]=Q1[i]+Q2[i]+Q3[i] 
如果投藥:假設投藥只准許於期初進行，令第 t 期期初所有將被投藥的網格集合為 V ，若於網格投
放撲殺率 X、範圍 R 的藥劑，將導致該網格所承接之前期期末病毒量 Q0 在該期不但無法成長,且無
法擴散給其相鄰網格，亦即 s=0;q=0
