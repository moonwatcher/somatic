#!/bin/bash

#somatic plot -p productive33+ -l 'c57bl6b1a'
#somatic plot -p productive33+ -l 'c57bl6fo'

somatic plot -p productive33+ -l 'c57bl6mzgf'
somatic plot -p productive33+ -l 'c57bl6mzspf'

somatic plot -p productive33+ -l 'c57bl6fogf'
somatic plot -p productive33+ -l 'c57bl6fospf'

somatic plot -p productive33+ -l 'c57bl6prebgf'
somatic plot -p productive33+ -l 'c57bl6prebspf'

somatic plot -p productive33+ -l 'c57bl6b1agf'
somatic plot -p productive33+ -l 'c57bl6b1aspf'

somatic plot -p productive33+ -l 'c57bl6gf'
somatic plot -p productive33+ -l 'c57bl6spf'

#c57bl6b1aspf 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3 vs c57bl6b1agf bc6fbf5a27a2812c6df549a1d9f2a8fb6b42fb6c
somatic plot-compare 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3 bc6fbf5a27a2812c6df549a1d9f2a8fb6b42fb6c

#c57bl6fospf 76319edd8ba6ed5a0363bdc10003304bac68da37 vs c57bl6fogf 37f4a13e5ea729a5b1ec11ecd8c8a5a552063f4c
somatic plot-compare 76319edd8ba6ed5a0363bdc10003304bac68da37 37f4a13e5ea729a5b1ec11ecd8c8a5a552063f4c

#c57bl6mzspf ab9a271b87361d47a6c52338d9731d903f8d83fc vs c57bl6mzgf 43863ec0ed3db08154a12118ed4fea862214ca28
somatic plot-compare ab9a271b87361d47a6c52338d9731d903f8d83fc 43863ec0ed3db08154a12118ed4fea862214ca28

#c57bl6prebspf d311f4e17bab01bc9aa083edee0e892f0d5b8765 vs c57bl6fospf 76319edd8ba6ed5a0363bdc10003304bac68da37
somatic plot-compare d311f4e17bab01bc9aa083edee0e892f0d5b8765 76319edd8ba6ed5a0363bdc10003304bac68da37

#c57bl6prebspf d311f4e17bab01bc9aa083edee0e892f0d5b8765 vs c57bl6mzspf ab9a271b87361d47a6c52338d9731d903f8d83fc
somatic plot-compare d311f4e17bab01bc9aa083edee0e892f0d5b8765 ab9a271b87361d47a6c52338d9731d903f8d83fc

#c57bl6prebspf d311f4e17bab01bc9aa083edee0e892f0d5b8765 vs c57bl6b1aspf 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3
somatic plot-compare d311f4e17bab01bc9aa083edee0e892f0d5b8765 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3

#c57bl6spf.675775d6886b7720f94083adc9f3564534f62c51 vs c57bl6gf.f8792812bb16acac61a313f04e83d3cdbc138336
somatic plot-compare 675775d6886b7720f94083adc9f3564534f62c51 f8792812bb16acac61a313f04e83d3cdbc138336

#c57bl6b1aspf.07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3 vs c57bl6b1agf.bc6fbf5a27a2812c6df549a1d9f2a8fb6b42fb6c
somatic plot-compare 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3 bc6fbf5a27a2812c6df549a1d9f2a8fb6b42fb6c

#c57bl6fospf.76319edd8ba6ed5a0363bdc10003304bac68da37 vs c57bl6b1aspf.07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3
somatic plot-compare 76319edd8ba6ed5a0363bdc10003304bac68da37 07b72597b7fba4b72825da3c0e2d9b5e1e70bfd3