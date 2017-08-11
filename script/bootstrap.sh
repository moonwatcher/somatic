#!/bin/bash
VERBOSITY="debug"
STRAIN="c57bl6"
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

rm -rf $BASE/html; mkdir $BASE/html
rm -rf $BASE/json; mkdir $BASE/json

mkdir -p $BASE/db/igblast/database
mkdir -p $BASE/db/igblast/optional_file

somatic rebuild --table gene

somatic -v $VERBOSITY genePopulate ../bootstrap/c57bl6/mouse_ighv_c57bl6.json
somatic -v $VERBOSITY geneAlign --strain "$STRAIN" --region VH --flanking 300 
somatic -v $VERBOSITY geneRss --strain "$STRAIN" --region VH --flanking 60 --distance 9
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region VH --title "IGHV $STRAIN" > $BASE/html/vh.html
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region VH --functionality F --title "IGHV $STRAIN Functional" > $BASE/html/vh_f.html
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region VH --functionality O --title "IGHV $STRAIN Open Reading Frame" > $BASE/html/vh_o.html
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region VH --functionality P --title "IGHV $STRAIN Pseudo" > $BASE/html/vh_p.html

somatic -v $VERBOSITY genePopulate ../bootstrap/c57bl6/mouse_ighd_c57bl6.json
somatic -v $VERBOSITY geneAlign --strain "$STRAIN" --region DH --flanking 300 
somatic -v $VERBOSITY geneRss --strain "$STRAIN" --region DH --flanking 60 --distance 9
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region DH --title "IGHD $STRAIN" > $BASE/html/dh.html

somatic -v $VERBOSITY genePopulate ../bootstrap/c57bl6/mouse_ighj_c57bl6.json
somatic -v $VERBOSITY geneAlign --strain "$STRAIN" --region JH --flanking 300 
somatic -v $VERBOSITY geneRss --strain "$STRAIN" --region JH --flanking 60 --distance 9
somatic -v $VERBOSITY geneHtml --strain "$STRAIN" --region JH --title "IGHJ $STRAIN" > $BASE/html/jh.html

somatic -v $VERBOSITY geneFasta --strain "$STRAIN" --region JH --profile aligned > $BASE/db/igblast/database/mouse_ighj
somatic -v $VERBOSITY geneFasta --strain "$STRAIN" --region VH --profile aligned > $BASE/db/igblast/database/mouse_ighv
somatic -v $VERBOSITY geneFasta --strain "$STRAIN" --region DH --profile aligned > $BASE/db/igblast/database/mouse_ighd

makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_jh -in $BASE/db/igblast/database/mouse_ighj
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_vh -in $BASE/db/igblast/database/mouse_ighv
makeblastdb -parse_seqids -dbtype nucl -input_type fasta -title mouse_dh -in $BASE/db/igblast/database/mouse_ighd

somatic -v $VERBOSITY geneIgblastAux --region VH --strain "$STRAIN" >  $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY geneIgblastAux --region DH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux
somatic -v $VERBOSITY geneIgblastAux --region JH --strain "$STRAIN" >> $BASE/db/igblast/optional_file/mouse_gl.aux

somatic -v $VERBOSITY geneInfo --region VH > $BASE/json/vh.json
somatic -v $VERBOSITY geneInfo --region DH > $BASE/json/dh.json 
somatic -v $VERBOSITY geneInfo --region JH > $BASE/json/jh.json 


# Library loading
somatic -v $VERBOSITY libraryPopulate ../bootstrap/c57bl6/library_c57bl6.json
