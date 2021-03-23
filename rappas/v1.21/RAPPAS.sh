#!/bin/bash

JAVA_OPTS="${JAVA_OPTS:=-Xmx2G}"
java $JAVA_OPTS -jar /usr/local/bin/RAPPAS.jar ${1+"$@"}
