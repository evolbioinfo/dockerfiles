#!/bin/bash

for REPO in pangolin pangoLEARN scorpio constellations pango-designation
do
    RELEASE=$(curl -s -H "Accept: application/vnd.github.v3+json"   https://api.github.com/repos/cov-lineages/$REPO/releases | jq '.[0].tag_name')
    echo $REPO ":" $RELEASE
done

for REPO in nextclade
do
    RELEASE=$(curl -s -H "Accept: application/vnd.github.v3+json"   https://api.github.com/repos/nextstrain/$REPO/releases | jq '.[0].tag_name')
    echo $REPO ":" $RELEASE
done
