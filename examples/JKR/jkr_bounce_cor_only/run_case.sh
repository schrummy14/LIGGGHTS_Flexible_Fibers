#!/bin/bash

rm -r post post_initial log* restart_*, v_out.txt;
../../../../src/lmp_auto -echo both < in.*

python post_process.py
