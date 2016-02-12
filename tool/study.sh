#!/bin/zsh

SESSION_NAME="somatic_study"
START_DIRECTORY="/Users/lg/thesis/study"

tmux new-session -c$START_DIRECTORY -s $SESSION_NAME -n "study" -d

tmux new-window -t"$SESSION_NAME:1" -n one "\
time somatic -v error study -p p33 -l c57bl6b01t01spfb1a --name c57bl6b01t01spfb1a_p33 --template VDJ > b01t01spfb1a_p33_vdj.csv; \
time somatic -v error study -p p33 -l c57bl6b02t01spfb1a --name c57bl6b02t01spfb1a_p33 --template VDJ > b02t01spfb1a_p33_vdj.csv; \
"

tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:2" -n two "\
time somatic -v error study -p p33 -l c57bl6b03t01spfb1a --name c57bl6b03t01spfb1a_p33 --template VDJ > b03t01spfb1a_p33_vdj.csv; \
time somatic -v error study -p p33 -l c57bl6b04t01spfb1a --name c57bl6b04t01spfb1a_p33 --template VDJ > b04t01spfb1a_p33_vdj.csv; \
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:3" -n three "\
time somatic -v error study -p p33 -l c57bl6b05t01spfb1a --name c57bl6b05t01spfb1a_p33 --template VDJ > b05t01spfb1a_p33_vdj.csv; \
time somatic -v error study -p p33 -l c57bl6b02t01gfb1a --name c57bl6b02t01gfb1a_p33 --template VDJ > b02t01gfb1a_p33_vdj.csv; \
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:4" -n four "\
time somatic -v error study -p p33 -l c57bl6b03t01gfb1a --name c57bl6b03t01gfb1a_p33 --template VDJ > b03t01gfb1a_p33_vdj.csv; \
time somatic -v error study -p p33 -l c57bl6b04t01gfb1a --name c57bl6b04t01gfb1a_p33 --template VDJ > b04t01gfb1a_p33_vdj.csv; \
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:5" -n five "\
time somatic -v error study -p p33 -l c57bl6b05t01gfb1a --name c57bl6b05t01gfb1a_p33 --template VDJ > b05t01gfb1a_p33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b01t01spfb1a --name c57bl6b01t01spfb1a_np33 --template VDJ > b01t01spfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b02t01spfb1a --name c57bl6b02t01spfb1a_np33 --template VDJ > b02t01spfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b03t01spfb1a --name c57bl6b03t01spfb1a_np33 --template VDJ > b03t01spfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b04t01spfb1a --name c57bl6b04t01spfb1a_np33 --template VDJ > b04t01spfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b05t01spfb1a --name c57bl6b05t01spfb1a_np33 --template VDJ > b05t01spfb1a_np33_vdj.csv; \
"
tmux new-window -c$START_DIRECTORY -t"$SESSION_NAME:6" -n six "\
time somatic -v error study -p np33 -l c57bl6b02t01gfb1a --name c57bl6b02t01gfb1a_np33 --template VDJ > b02t01gfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b03t01gfb1a --name c57bl6b03t01gfb1a_np33 --template VDJ > b03t01gfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b04t01gfb1a --name c57bl6b04t01gfb1a_np33 --template VDJ > b04t01gfb1a_np33_vdj.csv; \
time somatic -v error study -p np33 -l c57bl6b05t01gfb1a --name c57bl6b05t01gfb1a_np33 --template VDJ > b05t01gfb1a_np33_vdj.csv; \
"
tmux select-window -t "$SESSION_NAME:1"
tmux -2 attach-session -t "$SESSION_NAME"
