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
