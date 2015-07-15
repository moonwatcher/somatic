#!/bin/bash
VERBOSITY="debug"
STRAIN="C57BL/6"
BASE="/volume/albireo/canal/thesis"

rm -rf $BASE/html
rm -rf $BASE/db
cp -r ../db $BASE/

mkdir $BASE/html
mkdir -p $BASE/db/igblast/database
mkdir -p $BASE/db/igblast/optional_file

somatic rebuild --table gene

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighv_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r VH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r VH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --title "IGHV $STRAIN" > $BASE/html/vh.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality F --title "IGHV $STRAIN Functional" > $BASE/html/vh_f.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality O --title "IGHV $STRAIN Open Reading Frame" > $BASE/html/vh_o.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality P --title "IGHV $STRAIN Pseudo" > $BASE/html/vh_p.html

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighd_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r DH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r DH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r DH --title "IGHD $STRAIN" > $BASE/html/dh.html

somatic -v $VERBOSITY gene-populate ../bootstrap/c57bl6/mouse_ighj_c57bl6.json
somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r JH -F 300 
somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r JH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r JH --title "IGHJ $STRAIN" > $BASE/html/jh.html

somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r JH -p aligned > $BASE/db/igblast/database/mouse_c57bl6_ighj
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r VH -p aligned > $BASE/db/igblast/database/mouse_c57bl6_ighv
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r DH -p aligned > $BASE/db/igblast/database/mouse_c57bl6_ighd

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in $BASE/db/igblast/database/mouse_c57bl6_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in $BASE/db/igblast/database/mouse_c57bl6_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in $BASE/db/igblast/database/mouse_c57bl6_ighd

somatic -v $VERBOSITY gene-igblast-aux -r VH --strain "$STRAIN" >  $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r DH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r JH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux