# PASTdoc [PASTML Docker]

This docker file wraps [pastml](https://github.com/evolbioinfo/pastml) python3 module 
for Ancestor Character Reconstruction (ACR) and visualisation
on rooted phylogenetic trees.


Given a tree and its node annotations, it can either visualise them as-is, 
or infer ancestral node states based on the tip states. 

The states are visualised as different colours on the tree nodes using [Cytoscape.js](http://js.cytoscape.org/)


### How to run PASTdoc
As an input, one needs to provide a **rooted** phylogenetic tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
and a table containing tip states, 
in tab-delimited (by default) or csv format (to be specified with *--data_sep ,* option).

Basic usage:
```bash
docker run -v <path_to_the_folder_containing_the_tree_and_the_annotations>:/data:rw -t evolbioinfo/pastml --tree /data/<tree_file> --data /data/<annotation_file> --columns <one_or_more_column_names> --html_compressed /data/<map_name>
```


### Example
You can download [HIV1-A in Albania data](https://github.com/evolbioinfo/pastml/examples/Albania/data) as an example.
Let's assume that the tree and annotation files are in the Downloads folder, 
and are named respectively Albanian.tree.152tax.tre	and data.txt.

The data.txt is a comma-separated file, containing tip ids in the first column, 
and Country in the second column, i.e.:

id | Country
----- |  -----
98CMAJ6932 | Africa
98CMAJ6933 | Africa
96CMAJ6134 | Africa
00SEAY5240 | WestEurope
... | ...
02GRAY0303 | Greece
97YUAF9960 | EastEurope


To reconstruct and visualise the ancestral Country states, 
one needs to run the following command:

```bash
docker run -v ~/Downloads:/data:rw -t evolbioinfo/pastml --tree /data/Albanian.tree.152tax.tre --data /data/data.txt --data_sep , --columns Country --html_compressed /data/Albanian_map.html 
```

This will produce a file Albanian_map.html in the Downloads folder, 
that can be viewed with a browser.

###  Help

To see advanced options, run
```bash
docker run -t evolbioinfo/pastml -h
```
