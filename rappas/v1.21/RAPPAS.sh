#!/bin/bash
JAVA_OPTS= "${JAVA_OPTS:='-Xmx2g'}"

java $JAVA_OPTS -jar /usr/local/bin/RAPPAS/RAPPAS.jar ${1+"$@"}
