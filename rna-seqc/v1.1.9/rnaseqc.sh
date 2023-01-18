#!/usr/bin/env bash
java -Xmx6g -jar $( dirname -- "$0"; )/RNA-SeQC.jar ${1+"$@"}
