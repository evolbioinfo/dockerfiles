#!/usr/bin/env bash
java -Xmx6g -jar $( dirname -- "$0"; )/RNA-SeQC_v1.1.8.jar ${1+"$@"}
