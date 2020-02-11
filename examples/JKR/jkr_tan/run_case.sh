#!/bin/bash

rm -r post post_initial log* restart_*;
../../../../src/lmp_auto -echo both < in.*
