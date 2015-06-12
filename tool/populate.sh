#!/bin/bash
somatic rebuild --table gene
# somatic gene-populate ~/code/somatic/bootstrap/mouse_ighv_c57bl6_igblast.json
# somatic gene-align -r VH --strain 'C57BL/6' -F 300 
# somatic gene-rss -r VH --strain 'C57BL/6' -F 60 --distance 9
# somatic gene-html -r VH --strain 'C57BL/6' > ~/code/somatic/doc/vh_rss_igblast.html

somatic gene-populate ~/code/somatic/bootstrap/mouse_ighv_c57bl6_BN000872.json
somatic gene-align -r VH --strain 'C57BL/6' -F 300 
somatic gene-rss -r VH --strain 'C57BL/6' -F 60 --distance 9
somatic gene-html -r VH --strain 'C57BL/6' > ~/code/somatic/doc/vh_rss_BN000872.html

somatic gene-populate ~/code/somatic/bootstrap/mouse_ighd.json
somatic gene-align -r DH --strain 'C57BL/6' -F 300 
somatic gene-rss -r DH --strain 'C57BL/6' -F 60 --distance 9
somatic gene-html -r DH --strain 'C57BL/6' > ~/code/somatic/doc/dh_rss.html

somatic gene-populate ~/code/somatic/bootstrap/mouse_ighj.json
somatic gene-align -r JH --strain 'C57BL/6' -F 300 
somatic gene-rss -r JH --strain 'C57BL/6' -F 60 --distance 9
somatic gene-html -r JH --strain 'C57BL/6' > ~/code/somatic/doc/jh_rss.html
