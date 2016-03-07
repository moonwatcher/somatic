#!/bin/zsh

SESSION_NAME="somatic_study_vdj_fo"
START_DIRECTORY="/Users/lg/thesis/study_vdj_fo"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n "study_vdj_fo" -d

tmux new-window -t"$SESSION_NAME:1" -n one_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b02t02spffo'     --profile p33       --name 'spf_fo_b02t02_p33' > spf_fo_b02t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t01spffo'     --profile np33      --name 'spf_fo_b03t01_np33' > spf_fo_b03t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t01spffo'     --profile p33       --name 'spf_fo_b03t01_p33' > spf_fo_b03t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t02spffo'     --profile np33      --name 'spf_fo_b03t02_np33' > spf_fo_b03t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t02spffo'     --profile p33       --name 'spf_fo_b03t02_p33' > spf_fo_b03t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t03spffo'     --profile np33      --name 'spf_fo_b03t03_np33' > spf_fo_b03t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t03spffo'     --profile p33       --name 'spf_fo_b03t03_p33' > spf_fo_b03t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t01spffo'     --profile np33      --name 'spf_fo_b04t01_np33' > spf_fo_b04t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t01spffo'     --profile p33       --name 'spf_fo_b04t01_p33' > spf_fo_b04t01_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:2" -n two_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b04t02spffo'     --profile np33      --name 'spf_fo_b04t02_np33' > spf_fo_b04t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t02spffo'     --profile p33       --name 'spf_fo_b04t02_p33' > spf_fo_b04t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t01spffo'     --profile np33      --name 'spf_fo_b05t01_np33' > spf_fo_b05t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t01spffo'     --profile p33       --name 'spf_fo_b05t01_p33' > spf_fo_b05t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t02spffo'     --profile np33      --name 'spf_fo_b05t02_np33' > spf_fo_b05t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t02spffo'     --profile p33       --name 'spf_fo_b05t02_p33' > spf_fo_b05t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t03spffo'     --profile np33      --name 'spf_fo_b05t03_np33' > spf_fo_b05t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t03spffo'     --profile p33       --name 'spf_fo_b05t03_p33' > spf_fo_b05t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t04spffo'     --profile np33      --name 'spf_fo_b05t04_np33' > spf_fo_b05t04_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t04spffo'     --profile p33       --name 'spf_fo_b05t04_p33' > spf_fo_b05t04_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:3" -n three_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b02t01gffo'      --profile np33      --name 'gf_fo_b02t01_np33' > gf_fo_b02t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t01gffo'      --profile p33       --name 'gf_fo_b02t01_p33' > gf_fo_b02t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t02gffo'      --profile np33      --name 'gf_fo_b02t02_np33' > gf_fo_b02t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t02gffo'      --profile p33       --name 'gf_fo_b02t02_p33' > gf_fo_b02t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t03gffo'      --profile np33      --name 'gf_fo_b02t03_np33' > gf_fo_b02t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t03gffo'      --profile p33       --name 'gf_fo_b02t03_p33' > gf_fo_b02t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t01gffo'      --profile np33      --name 'gf_fo_b03t01_np33' > gf_fo_b03t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t01gffo'      --profile p33       --name 'gf_fo_b03t01_p33' > gf_fo_b03t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t02gffo'      --profile np33      --name 'gf_fo_b03t02_np33' > gf_fo_b03t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t02gffo'      --profile p33       --name 'gf_fo_b03t02_p33' > gf_fo_b03t02_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:4" -n four_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b03t03gffo'      --profile np33      --name 'gf_fo_b03t03_np33' > gf_fo_b03t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t03gffo'      --profile p33       --name 'gf_fo_b03t03_p33' > gf_fo_b03t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t04gffo'      --profile np33      --name 'gf_fo_b03t04_np33' > gf_fo_b03t04_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b03t04gffo'      --profile p33       --name 'gf_fo_b03t04_p33' > gf_fo_b03t04_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t01gffo'      --profile np33      --name 'gf_fo_b04t01_np33' > gf_fo_b04t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t01gffo'      --profile p33       --name 'gf_fo_b04t01_p33' > gf_fo_b04t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t02gffo'      --profile np33      --name 'gf_fo_b04t02_np33' > gf_fo_b04t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t02gffo'      --profile p33       --name 'gf_fo_b04t02_p33' > gf_fo_b04t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t03gffo'      --profile np33      --name 'gf_fo_b04t03_np33' > gf_fo_b04t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t03gffo'      --profile p33       --name 'gf_fo_b04t03_p33' > gf_fo_b04t03_p33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:5" -n five_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b04t04gffo'      --profile np33      --name 'gf_fo_b04t04_np33' > gf_fo_b04t04_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b04t04gffo'      --profile p33       --name 'gf_fo_b04t04_p33' > gf_fo_b04t04_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t01gffo'      --profile np33      --name 'gf_fo_b05t01_np33' > gf_fo_b05t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t01gffo'      --profile p33       --name 'gf_fo_b05t01_p33' > gf_fo_b05t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t02gffo'      --profile np33      --name 'gf_fo_b05t02_np33' > gf_fo_b05t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t02gffo'      --profile p33       --name 'gf_fo_b05t02_p33' > gf_fo_b05t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t03gffo'      --profile np33      --name 'gf_fo_b05t03_np33' > gf_fo_b05t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b05t03gffo'      --profile p33       --name 'gf_fo_b05t03_p33' > gf_fo_b05t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b01t01spffo'     --profile np33      --name 'spf_fo_b01t01_np33' > spf_fo_b01t01_np33.csv;\
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:6" -n six_vdj_fo "\
somatic study --preset VDJ --library 'c57bl6b01t01spffo'     --profile p33       --name 'spf_fo_b01t01_p33' > spf_fo_b01t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b01t02spffo'     --profile np33      --name 'spf_fo_b01t02_np33' > spf_fo_b01t02_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b01t02spffo'     --profile p33       --name 'spf_fo_b01t02_p33' > spf_fo_b01t02_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b01t03spffo'     --profile np33      --name 'spf_fo_b01t03_np33' > spf_fo_b01t03_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b01t03spffo'     --profile p33       --name 'spf_fo_b01t03_p33' > spf_fo_b01t03_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t01spffo'     --profile np33      --name 'spf_fo_b02t01_np33' > spf_fo_b02t01_np33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t01spffo'     --profile p33       --name 'spf_fo_b02t01_p33' > spf_fo_b02t01_p33.csv;\
somatic study --preset VDJ --library 'c57bl6b02t02spffo'     --profile np33      --name 'spf_fo_b02t02_np33' > spf_fo_b02t02_np33.csv;\
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"
