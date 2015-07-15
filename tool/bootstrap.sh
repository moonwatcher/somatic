#!/bin/bash
VERBOSITY="debug"
STRAIN="C57BL/6"
BASE="/volume/albireo/canal/thesis"

[ -d "$BASE/db" ] || mkdir -p "$BASE/db"
[ -d "$BASE/db/igblast" ] && rm -rf "$BASE/db/igblast"
cp -r ../db/igblast "$BASE/db/"

if [ ! -f "$BASE/db/mus_musculus.grcm38.dna.chromosome.12.fa" ]; then
    (
        cd $BASE/db
        curl -O "http://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.12.fa.gz"
        gzip -d -c Mus_musculus.GRCm38.dna.chromosome.12.fa.gz > mus_musculus.grcm38.dna.chromosome.12.fa
        rm -f Mus_musculus.GRCm38.dna.chromosome.12.fa.gz
    )
fi

rm -rf $BASE/html
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

somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r JH -p aligned > $BASE/db/igblast/database/mouse_ighj
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r VH -p aligned > $BASE/db/igblast/database/mouse_ighv
somatic -v $VERBOSITY gene-fasta --strain "$STRAIN" -r DH -p aligned > $BASE/db/igblast/database/mouse_ighd

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_jh -in $BASE/db/igblast/database/mouse_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_vh -in $BASE/db/igblast/database/mouse_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_imgt_dh -in $BASE/db/igblast/database/mouse_ighd

somatic -v $VERBOSITY gene-igblast-aux -r VH --strain "$STRAIN" >  $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r DH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY gene-igblast-aux -r JH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux