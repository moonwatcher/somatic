#!/bin/bash
VERBOSITY="debug"
STRAIN="C57BL/6"

rm -rf ../html
rm -rf ../db/igblast/database
rm -rf ../db/igblast/optional_file

mkdir ../html
mkdir ../db/igblast/database
mkdir ../db/igblast/optional_file

somatic rebuild --table gene

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighv_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r VH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r VH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --title "IGHV $STRAIN" > ../html/vh.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality F --title "IGHV $STRAIN Functional" > ../html/vh_f.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality O --title "IGHV $STRAIN Open Reading Frame" > ../html/vh_o.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality P --title "IGHV $STRAIN Pseudo" > ../html/vh_p.html

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighd_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r DH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r DH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r DH --title "IGHD $STRAIN" > ../html/dh.html

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighj_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r JH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r JH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r JH --title "IGHJ $STRAIN" > ../html/jh.html

somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r JH -p aligned > ../db/igblast/database/mouse_c57bl6_ighj
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r VH -p aligned > ../db/igblast/database/mouse_c57bl6_ighv
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r DH -p aligned > ../db/igblast/database/mouse_c57bl6_ighd

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in ../db/igblast/database/mouse_c57bl6_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in ../db/igblast/database/mouse_c57bl6_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in ../db/igblast/database/mouse_c57bl6_ighd

somatic -v $VERBOSITY gene-igblast-aux -r VH --strain "$STRAIN" >  ../db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r DH --strain "$STRAIN" >> ../db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r JH --strain "$STRAIN" >> ../db/igblast/optional_file/mouse_gl.aux