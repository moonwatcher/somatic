#!/bin/zsh

SESSION_NAME="somatic_study_vdj_b1a"
START_DIRECTORY="/Users/lg/thesis/study_vdj_b1a"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n "study_vdj_b1a" -d

tmux new-window -t"$SESSION_NAME:1" -n one_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b02t01gfb1a'      --profile np33      --name 'gf_b1a_b02_np33' > gf_b1a_b02_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b02t01gfb1a'      --profile p33       --name 'gf_b1a_b02_p33' > gf_b1a_b02_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b03t01gfb1a'      --profile np33      --name 'gf_b1a_b03_np33' > gf_b1a_b03_np33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:2" -n two_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b03t01gfb1a'      --profile p33       --name 'gf_b1a_b03_p33' > gf_b1a_b03_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b04t01gfb1a'      --profile np33      --name 'gf_b1a_b04_np33' > gf_b1a_b04_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b04t01gfb1a'      --profile p33       --name 'gf_b1a_b04_p33' > gf_b1a_b04_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:3" -n three_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b05t01gfb1a'      --profile np33      --name 'gf_b1a_b05_np33' > gf_b1a_b05_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b05t01gfb1a'      --profile p33       --name 'gf_b1a_b05_p33' > gf_b1a_b05_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b01t01spfb1a'     --profile np33      --name 'spf_b1a_b01_np33' > spf_b1a_b01_np33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:4" -n four_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b01t01spfb1a'     --profile p33       --name 'spf_b1a_b01_p33' > spf_b1a_b01_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b02t01spfb1a'     --profile np33      --name 'spf_b1a_b02_np33' > spf_b1a_b02_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b02t01spfb1a'     --profile p33       --name 'spf_b1a_b02_p33' > spf_b1a_b02_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:5" -n five_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b03t01spfb1a'     --profile np33      --name 'spf_b1a_b03_np33' > spf_b1a_b03_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b03t01spfb1a'     --profile p33       --name 'spf_b1a_b03_p33' > spf_b1a_b03_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b04t01spfb1a'     --profile np33      --name 'spf_b1a_b04_np33' > spf_b1a_b04_np33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:6" -n six_vdj_b1a "\
somatic -v error study --preset VDJ --library 'c57bl6b04t01spfb1a'     --profile p33       --name 'spf_b1a_b04_p33' > spf_b1a_b04_p33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b05t01spfb1a'     --profile np33      --name 'spf_b1a_b05_np33' > spf_b1a_b05_np33.csv;\
somatic -v error study --preset VDJ --library 'c57bl6b05t01spfb1a'     --profile p33       --name 'spf_b1a_b05_p33' > spf_b1a_b05_p33.csv;\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"


