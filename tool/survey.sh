#!/bin/zsh

SESSION_NAME="somatic_survey"
START_DIRECTORY="/Users/lg/thesis/plot"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n etc -d
new-window -t"$SESSION_NAME:1" -n one "\
somatic -v debug survey --library 'c57bl6b1agf'      --profile np33      --name 'gf_b1a_np33';\
somatic -v debug survey --library 'c57bl6b1agf'      --profile p33       --name 'gf_b1a_p33';\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile np33      --name 'spf_b1a_np33';\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile npn033    --name 'spf_b1a_npn033';\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile npn133    --name 'spf_b1a_npn133';\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile p33       --name 'spf_b1a_p33';\
"
tmux new-window -t"$SESSION_NAME:2" -n two "\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile pn033     --name 'spf_b1a_pn033';\
somatic -v debug survey --library 'c57bl6b1aspf'     --profile pn133     --name 'spf_b1a_pn133';\
somatic -v debug survey --library 'c57bl6fogf'       --profile np33      --name 'gf_fo_np33';\
somatic -v debug survey --library 'c57bl6fogf'       --profile p33       --name 'gf_fo_p33';\
somatic -v debug survey --library 'c57bl6fospf'      --profile np33      --name 'spf_fo_np33';\
somatic -v debug survey --library 'c57bl6fospf'      --profile p33       --name 'spf_fo_p33';\
"
tmux new-window -t"$SESSION_NAME:3" -n three "\
somatic -v debug survey --library 'c57bl6gf'         --profile np33      --name 'gf_np33';\
somatic -v debug survey --library 'c57bl6gf'         --profile p33       --name 'gf_p33';\
somatic -v debug survey --library 'c57bl6mzgf'       --profile np33      --name 'gf_mz_np33';\
somatic -v debug survey --library 'c57bl6mzgf'       --profile p33       --name 'gf_mz_p33';\
somatic -v debug survey --library 'c57bl6mzspf'      --profile np33      --name 'spf_mz_np33';\
somatic -v debug survey --library 'c57bl6mzspf'      --profile p33       --name 'spf_mz_p33';\
"
tmux new-window -t"$SESSION_NAME:4" -n four "\
somatic -v debug survey --library 'c57bl6prebgf'     --profile np33      --name 'gf_preb_np33';\
somatic -v debug survey --library 'c57bl6prebgf'     --profile p33       --name 'gf_preb_p33';\
somatic -v debug survey --library 'c57bl6prebspf'    --profile np33      --name 'spf_preb_np33';\
somatic -v debug survey --library 'c57bl6prebspf'    --profile p33       --name 'spf_preb_p33';\
somatic -v debug survey --library 'c57bl6spf'        --profile np33      --name 'spf_np33';\
somatic -v debug survey --library 'c57bl6spf'        --profile p33       --name 'spf_p33';\
"
tmux new-window -t"$SESSION_NAME:5" -n five "\
somatic -v debug survey --library 'c57bl6spfiabm'    --profile np33      --name 'spf_iabm_np33';\
somatic -v debug survey --library 'c57bl6spfiabm'    --profile p33       --name 'spf_iabm_p33';\
somatic -v debug survey --library 'c57bl6spfiaspl'   --profile np33      --name 'spf_iaspl_np33';\
somatic -v debug survey --library 'c57bl6spfiaspl'   --profile p33       --name 'spf_iaspl_p33';\
somatic -v debug survey --library 'c57bl6spfiftl'    --profile np33      --name 'spf_iftl_np33';\
"
tmux new-window -t"$SESSION_NAME:6" -n six "\
somatic -v debug survey --library 'c57bl6spfiftl'    --profile p33       --name 'spf_iftl_p33';\
somatic -v debug survey --library 'c57bl6spfprebabm' --profile np33      --name 'spf_prebabm_np33';\
somatic -v debug survey --library 'c57bl6spfprebabm' --profile p33       --name 'spf_prebabm_p33';\
somatic -v debug survey --library 'c57bl6spfprebftl' --profile np33      --name 'spf_prebftl_np33';\
somatic -v debug survey --library 'c57bl6spfprebftl' --profile p33       --name 'spf_prebftl_p33';\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"
