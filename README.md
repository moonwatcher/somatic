Generate an igblast database
============================

```
somatic load-reference imgt/mouse/mouse_igmt_dh.fasta imgt/mouse/mouse_igmt_jh.fasta imgt/mouse/mouse_igmt_vh.fasta
somatic to-blast-fasta -r V > igblast/database/mouse_imgt_vh
somatic to-blast-fasta -r D > igblast/database/mouse_imgt_dh
somatic to-blast-fasta -r J > igblast/database/mouse_imgt_jh

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in igblast/database/mouse_imgt_vh
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in igblast/database/mouse_imgt_dh
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in igblast/database/mouse_imgt_jh

somatic to-auxiliary -r V >> igblast/optional_file/mouse_gl.aux
somatic to-auxiliary -r D >> igblast/optional_file/mouse_gl.aux
somatic to-auxiliary -r J >> igblast/optional_file/mouse_gl.aux
```

Load a sample
=============
```gzcat A4N54_l01n01_1FObB61.fastq.gz|somatic populate -l FObB61```