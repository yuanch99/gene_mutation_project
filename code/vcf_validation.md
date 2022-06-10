### D1243

Tissue Format: Frozen

**Old file:**

**TMB: 496.24**

D1243.hc.hrdflt.vcf_ext.tsv.annovar.out_FINAL_rev26.5.tsv
PASS: 105202
mutect2.D1243.MMR100B1.b37.merged.vcf 
PASS: 25314
D1243.hc.hrdflt.vcf
PASS: 94861

<u>POLE:</u>
p.S461P c.T1381C
12	133249842	133249842	A	G 0.353

p.I1423I c.C4269T
12      133220444       133220444       G       A 0.087

<u>POLD1:</u>
c.C111T:p.D37D
c.C1449T:p.G483G
c.C2106T:p.G702G

**New file:**

**TMB: 42.49**

D1243__MMR100B1.mutect2.filtered.wes.vcf.gz
PASS: 2053

<u>POLE:</u>
p.S461
12  132673256 132673256 A	G 0.347
*orientation*

p.I1423I c.C4269T
chr12   132643858       .       G       A       .       PASS

<u>POLD1:</u>
p.D37D 19:50398962
*orientation*
p.G483G 19:50406472
*orientation*
p.G702G 19:50409618
*orientation*

<img src="/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/vcf_validation/output/D1243.png" alt="D1243" style="zoom:24%;" />


```bash
(base) yuanchang@YuandeMacBook-Air data % cat mutect2.D1243.MMR100B1.b37.merged.vcf | grep -v '^#' | grep 'PASS' | wc -l
   26243
(base) yuanchang@YuandeMacBook-Air data % gzcat D1243__MMR100B1.mutect2.filtered.wes.vcf.gz | grep -v '^#' | grep 'orientation\|PASS' | wc -l
   26302
(base) yuanchang@YuandeMacBook-Air data % gzcat D1243__MMR100B1.mutect2.filtered.wes.vcf.gz | grep -v '^#' | grep 'PASS' | wc -l 
    2629
```



### MD1341T3

Tissue Format: Frozen

**Old file:**

**TMB: 282**

MD1341T3.hc.hrdflt.vcf
PASS: 140457
MD1341T3.hc.hrdflt.vcf_fixed.tsv.annovar.out_FINAL_rev26.2.tsv_withTargetInfo.tsv
PASS: 151668
mutect2.MD1341T3.MD1341B1.with_PON.merged.20171110.vcf
PASS: 30050


<u>POLE:</u>

c.C4098A:p.F1366L
12   133225566

c.A2590G:p.I864V
12   133240706

c.C1690A:p.P564T
12   133248905

==c890T:p.S297F==
12   133253151

<u>POLD1:</u>



**New file:**

**TMB:17.07**

MD1341T3__MD1341B1.mutect2.all.Somatic.annotated-snpeff.wes.vcf.gz
PASS: 1225
+orientation: 25435

<u>POLE</u>

G T:p.F1366L
12 132648980
*orientation*

T C: p.I864V
12 132664120
*orientation*

c.C1690A:p.P564T
None

==G A: p.S297F==
12 132676565
*orientation* 0.43

<u>POLD1</u>

chr12  12726272    .    T    TA



### JG



**Old file**

**TMBL: 181.66**

JG.hc.hrdflt.vcf_ext.tsv.annovar.out_FINAL_rev26.5.tsv
PASS: 135176
mutect2.JG.MMR152B1.b37.merged.vcf
PASS: 16060
JG.hc.hrdflt.vcf
PASS: 121706

<u>POLE</u>

c.A6252G:p.S2084S

c.A4530G:p.A1510A

c.G3156A:p.T1052T

<u>POLD1</u>

c.G952A:p.E318K

19 50905980    G    A



**New file:**

**TMB:179.9**

JG__MMR152B1.mutect2.filtered.wes.vcf.gz
PASS: 9603
+orientation: 11519
New 

<u>POLE</u>
c.357T>G|p.Phe119Leu
chr12   132680020       .       A       C  

<u>POLD1</u>

c.952G>A|p.Glu318Lys

chr19  50402723    .    G    A



### Notes

37 -> 38

Seq Who Where 

- sample 
  - which are cor and which are not (why) (format age-sample age-patient new/old)
  - TMB diff
  - Old vcf source (tmb)
  - seq source



read pile up IGV of 1



para vs seq



H



#files

'Good files' : <10%

"MD1460T1","MD1440T2","MD1515T1","MD1564T1","4","6","MD1456T2","MD1630T05","MD1468T1","MD1630T02","JG"

'Bad files' : >10%
"D1807_1","D1807_2","MD1468T2","D1577","MD1566T2","D132","MD1385T2_2","D1806","MD1454T1","MD1456T3","D1120","D1804_3","D1144","MD1475T1","MD1385T2_1","MD0134T10","D1121","X3","D1243","D1804_2","D1423","X2","MD1293T1","D1410","MD1273T25","MD0134T11","MD190T2","MD101T10","MD1341T3","MD1273T26","D1764","M1273T2","MD0134T12","MD1303T1"

==MD1475T1__MD1475B1== 811310 snps???? wtf



clear signature (e.g. C to A SNP that occurs predominantly for the middle C in the DNA sequence CCG)



