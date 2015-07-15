#!/bin/bash

somatic survey -p productive33+ -l 'c57bl6fospf' --name 'spffo33'
somatic survey -p productive33+ -l 'c57bl6mzspf' --name 'spfmz33'
somatic survey -p productive33+ -l 'c57bl6prebspf' --name 'spfpreb33'
somatic survey -p productive33+ -l 'c57bl6b1aspf' --name 'spfb1a33'

somatic survey -p productive33+ -l 'c57bl6fogf' --name 'gffo33'
somatic survey -p productive33+ -l 'c57bl6mzgf' --name 'gfmz33'
somatic survey -p productive33+ -l 'c57bl6prebgf' --name 'gfpreb33'
somatic survey -p productive33+ -l 'c57bl6b1agf' --name 'gfb1a33'

# somatic survey -p productive33+ -l 'c57bl6gf' --name 'gf33'
# somatic survey -p productive33+ -l 'c57bl6spf' --name 'spf33'

somatic survey -p productive33+ -l 'c57bl6spfiabm33' --name 'spfiabm33'
somatic survey -p productive33+ -l 'c57bl6spfiaspl33' --name 'spfiaspl33'
somatic survey -p productive33+ -l 'c57bl6spfiftl33' --name 'spfiftl33'
somatic survey -p productive33+ -l 'c57bl6spfprebabm33' --name 'spfprebabm33'
somatic survey -p productive33+ -l 'c57bl6spfprebftl33' --name 'spfprebftl33'

somatic survey -p productiveN033+ -l 'c57bl6b1aspf' --name 'spfb1a33n0'
somatic survey -p productiveN1+33+ -l 'c57bl6b1aspf' --name 'spfb1a33n1+'

# single
somatic plot spfiabm33
somatic plot spfiaspl33
somatic plot spfiftl33
somatic plot spfprebabm33
somatic plot spfprebftl33

# comparisons
somatic plot spfb1a33 gfb1a33
somatic plot spffo33 gffo33
somatic plot spfmz33 gfmz33
somatic plot spfpreb33 spffo33
somatic plot spfpreb33 spfmz33
somatic plot spfpreb33 spfb1a33
somatic plot spf33 gf33
somatic plot spfb1a33 gfb1a33
somatic plot spffo33 spfb1a33

somatic plot spffo33 spfb1a33n0
somatic plot spfb1a33n1+ spfb1a33n0

somatic plot spffo33 spfiftl33
somatic plot spfb1a33 spfiftl33
somatic plot spfb1a33n0 spfiftl33
somatic plot spfb1a33n1+ spfiftl33