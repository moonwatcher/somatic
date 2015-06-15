#!/bin/bash

somatic survey -p productive33+ -l 'c57bl6fospf' --name 'spffo33'
somatic survey -p productive33+ -l 'c57bl6mzspf' --name 'spfmz33'
somatic survey -p productive33+ -l 'c57bl6prebspf' --name 'spfpreb33'
somatic survey -p productive33+ -l 'c57bl6b1aspf' --name 'spfb1a33'

somatic survey -p productive33+ -l 'c57bl6fogf' --name 'gffo33'
somatic survey -p productive33+ -l 'c57bl6mzgf' --name 'gfmz33'
somatic survey -p productive33+ -l 'c57bl6prebgf' --name 'gfpreb33'
somatic survey -p productive33+ -l 'c57bl6b1agf' --name 'gfb1a33'

somatic survey -p productive33+ -l 'c57bl6gf' --name 'gf33'
somatic survey -p productive33+ -l 'c57bl6spf' --name 'spf33'

somatic plot spfb1a33 gfb1a33
somatic plot spffo33 gffo33
somatic plot spfmz33 gfmz33
somatic plot spfpreb33 spffo33
somatic plot spfpreb33 spfmz33
somatic plot spfpreb33 spfb1a33
somatic plot spf33 gf33
somatic plot spfb1a33 gfb1a33
somatic plot spffo33 spfb1a33
