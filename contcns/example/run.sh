#!/bin/bash

oc2mkdb wrk_dir ont_read_list.txt
oc2pm -t 4 wrk_dir candidates.txt
