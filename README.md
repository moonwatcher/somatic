IMGT resources
==============
```
http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+IGHV&species=Mus+musculus
http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+IGHD&species=Mus+musculus
http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+IGHJ&species=Mus+musculus
```


Generate an igblast database
============================

```
somatic load-reference imgt/mouse/mouse_imgt_dh.fasta imgt/mouse/mouse_imgt_jh.fasta imgt/mouse/mouse_imgt_vh.fasta
somatic to-blast-fasta -r VH > igblast/database/mouse_imgt_vh
somatic to-blast-fasta -r DH > igblast/database/mouse_imgt_dh
somatic to-blast-fasta -r JH > igblast/database/mouse_imgt_jh

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in igblast/database/mouse_imgt_vh
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in igblast/database/mouse_imgt_dh
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in igblast/database/mouse_imgt_jh

somatic to-auxiliary -r VH >> igblast/optional_file/mouse_gl.aux
somatic to-auxiliary -r DH >> igblast/optional_file/mouse_gl.aux
somatic to-auxiliary -r JH>> igblast/optional_file/mouse_gl.aux
```

Load a sample
=============
```gzcat A4N54_l01n01_1FObB61.fastq.gz|somatic populate -l FObB61```
