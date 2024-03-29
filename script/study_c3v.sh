#!/bin/zsh

SESSION_NAME="somatic_study_c3v"
START_DIRECTORY="/Users/lg/thesis/study_c3v"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n "study_c3v" -d

tmux new-window -t"$SESSION_NAME:1" -n one_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6b1agf'      --profile np33      --name 'gf_b1a_np33_c3v' > gf_b1a_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1agf'      --profile p33       --name 'gf_b1a_p33_c3v' > gf_b1a_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile np33      --name 'spf_b1a_np33_c3v' > spf_b1a_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile npn033    --name 'spf_b1a_npn033_c3v' > spf_b1a_npn033_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile npn133    --name 'spf_b1a_npn133_c3v' > spf_b1a_npn133_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile p33       --name 'spf_b1a_p33_c3v' > spf_b1a_p33_c3v.csv;\
"

tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:2" -n two_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile pn033     --name 'spf_b1a_pn033_c3v' > spf_b1a_pn033_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6b1aspf'     --profile pn133     --name 'spf_b1a_pn133_c3v' > spf_b1a_pn133_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6fogf'       --profile np33      --name 'gf_fo_np33_c3v' > gf_fo_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6fogf'       --profile p33       --name 'gf_fo_p33_c3v' > gf_fo_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6fospf'      --profile np33      --name 'spf_fo_np33_c3v' > spf_fo_np33_c3v.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:3" -n three_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6fospf'      --profile p33       --name 'spf_fo_p33_c3v' > spf_fo_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6mzgf'       --profile np33      --name 'gf_mz_np33_c3v' > gf_mz_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6mzgf'       --profile p33       --name 'gf_mz_p33_c3v' > gf_mz_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6mzspf'      --profile np33      --name 'spf_mz_np33_c3v' > spf_mz_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6mzspf'      --profile p33       --name 'spf_mz_p33_c3v' > spf_mz_p33_c3v.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:4" -n four_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6prebgf'     --profile np33      --name 'gf_preb_np33_c3v' > gf_preb_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6prebgf'     --profile p33       --name 'gf_preb_p33_c3v' > gf_preb_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6prebspf'    --profile np33      --name 'spf_preb_np33_c3v' > spf_preb_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6prebspf'    --profile p33       --name 'spf_preb_p33_c3v' > spf_preb_p33_c3v.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:5" -n five_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6spfiabm'    --profile np33      --name 'spf_iabm_np33_c3v' > spf_iabm_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfiabm'    --profile p33       --name 'spf_iabm_p33_c3v' > spf_iabm_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfiaspl'   --profile np33      --name 'spf_iaspl_np33_c3v' > spf_iaspl_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfiaspl'   --profile p33       --name 'spf_iaspl_p33_c3v' > spf_iaspl_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfiftl'    --profile np33      --name 'spf_iftl_np33_c3v' > spf_iftl_np33_c3v.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:6" -n six_c3v "\
time somatic -v error study --preset C3V --library 'c57bl6spfiftl'    --profile p33       --name 'spf_iftl_p33_c3v' > spf_iftl_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfprebabm' --profile np33      --name 'spf_prebabm_np33_c3v' > spf_prebabm_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfprebabm' --profile p33       --name 'spf_prebabm_p33_c3v' > spf_prebabm_p33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfprebftl' --profile np33      --name 'spf_prebftl_np33_c3v' > spf_prebftl_np33_c3v.csv;\
time somatic -v error study --preset C3V --library 'c57bl6spfprebftl' --profile p33       --name 'spf_prebftl_p33_c3v' > spf_prebftl_p33_c3v.csv;\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"


