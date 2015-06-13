#!/bin/bash
VERBOSITY="debug"
STRAIN="C57BL/6"

# somatic rebuild --table gene
# somatic -v $VERBOSITY gene-populate ~/code/somatic/bootstrap/mouse_ighv_c57bl6_igblast.json
# somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r VH -F 300 
# somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r VH -F 60 --distance 9
# somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH > ~/code/somatic/doc/vh_rss_igblast.html

# somatic -v $VERBOSITY gene-populate ~/code/somatic/bootstrap/mouse_ighv_c57bl6_bn000872.json
# somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r VH -F 300 
# somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r VH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --title "IGHV $STRAIN BN000872" > ~/code/somatic/doc/vh_rss_bn000872.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality F --title "IGHV C57BL/6 Functional BN000872" > ~/code/somatic/doc/vh_rss_bn000872_f.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality O --title "IGHV C57BL/6 Open Reading Frame BN000872" > ~/code/somatic/doc/vh_rss_bn000872_o.html
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r VH --functionality P --title "IGHV C57BL/6 Pseudo BN000872" > ~/code/somatic/doc/vh_rss_bn000872_p.html

# somatic -v $VERBOSITY gene-populate ~/code/somatic/bootstrap/mouse_ighd.json
# somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r DH -F 300 
# somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r DH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r DH --title "IGHD $STRAIN" > ~/code/somatic/doc/dh_rss.html

# somatic -v $VERBOSITY gene-populate ~/code/somatic/bootstrap/mouse_ighj.json
# somatic -v $VERBOSITY gene-align --strain "$STRAIN" -r JH -F 300 
# somatic -v $VERBOSITY gene-rss --strain "$STRAIN" -r JH -F 60 --distance 9
somatic -v $VERBOSITY gene-html --strain "$STRAIN" -r JH --title "IGHJ $STRAIN" > ~/code/somatic/doc/jh_rss.html
