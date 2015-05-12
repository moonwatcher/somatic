Bootstrap the system
====================

**Build the mongodb indexes**
```
somatic rebuild
```

**Load genes for C57BL/6**
```
somatic gene-populate bootstrap/mouse_ighd.json 
somatic gene-populate bootstrap/mouse_ighj.json 
somatic gene-populate bootstrap/mouse_ighv_c57bl6_bn000872.json 
somatic gene-populate bootstrap/mouse_ighv_extra.json 
```

**Align the genes to the reference with up to 200bp flanking sequences**
```
somatic gene-align -F 200 --strain 'C57BL/6'
```

**Create FASTA files with the C57BL/6 genes**
```
somatic gene-fasta --strain 'C57BL/6' -r JH -p aligned > db/igblast/database/mouse_c57bl6_ighj
somatic gene-fasta --strain 'C57BL/6' -r VH -p aligned > db/igblast/database/mouse_c57bl6_ighv
somatic gene-fasta --strain 'C57BL/6' -r DH -p aligned > db/igblast/database/mouse_c57bl6_ighd
```

**Generate an igblast database**
```
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in db/igblast/database/mouse_c57bl6_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in db/igblast/database/mouse_c57bl6_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in db/igblast/database/mouse_c57bl6_ighd

somatic gene-igblast-aux -r VH --strain 'C57BL/6' >  db/igblast/optional_file/mouse_gl.aux
somatic gene-igblast-aux -r DH --strain 'C57BL/6' >> db/igblast/optional_file/mouse_gl.aux
somatic gene-igblast-aux -r JH --strain 'C57BL/6' >> db/igblast/optional_file/mouse_gl.aux
```

Load a sample
=============
```zcat "A4G6U l01n01 b6_spf_preb_2.fastq.gz" | somatic populate --strain C57BL/6 --library c57bl6b02t01spfpreb```

Look at an alignment diagram
============================
```
somatic view -L 2 -S 128 -p productive

R    name        F  c57bl6b03t03spffo | @M00595:21:000000000-A4G6U:1:1101:10492:4386 1:N:0:ATCAC | VH : J558.16.106 | DH : DSP2.2 | JH : IGHJ1 | productive | Q 36.25
                    0  2   5   8   11  14  17  20  23  26  29  32  35  38  41  44  47  50  53  56  59  62  65  68  71  74  77  80  83  86  89  92  95  98  101 104 107 110 113 116 119 122 125 128 131 134 137 140 143 
                    CC AGC ACA GCC TAC ATG GAG CTC CGC AGC CTG ACA TCT GAG GAC ACT GCA GTC TAT TAC TGT GCA AGA CGG GGT TAT GAT TAC GAC GGT TAC TGG TAC TTC GAT GTC TGG GGC ACA GGG ACC ACG GTC ACC GTC TCC TCA GGT AAG 
                    CD DFF HG< F22 F>@ 2HH G0C B>C DGF HHH HGH HGD B3F G4? F?? BCH HHH GFH HHF HHH HGE HHG GE> GGH HHH HHH HFF F>E EEA FFG HHH HHH HFG HFF HHH GGA GDH HFF GFG HHH HGH GGG HGH GGG GGG GGG GEF FCD FFF 
                       S   T   A   Y   M   E   L   R   S   L   T   S   E   D   T   A   V   Y   Y   C   A   R   R   G   Y   D   Y   D   G   Y   W   Y   F   D   V   W   G   T   G   T   T   V   T   V   S   S   G   K   
VH   J558.16.106 F  -- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
DH   DSP2.2      F                                                                                                     --- --- --- --- 
JH   IGHJ1       F                                                                                                                         --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- - 
V-D                                                                                                            --- --- 
CDR3                                                                                               --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
D-J                                                                                                                                    --- 
                                                                                                                                       PNN 
```
