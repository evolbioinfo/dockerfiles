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
Let's assume that the tree and annotation files are in the Downloads folder, 
and are named respectively tree.nwk and states.csv.

The states.csv is a comma-separated file, containing tip ids in the first column, 
and several named columns, including *Location*, i.e.:


Tip_id | ... | Location | ...
----- |  ----- | ----- | -----
1 | ... | Africa | ...
2 | ... | Asia | ...
3 | ... | Africa | ...
... | ... | ... | ...


To reconstruct and visualise the ancestral Location states, 
one needs to run the following command:

```bash
docker run -v ~/Downloads:/data:rw -t evolbioinfo/pastml --tree /data/tree.nwk --data /data/states.csv --data_sep , --columns Location --html_compressed /data/location_map.html
```

This will produce a file location_map.html in the Downloads folder, 
that can be viewed with a browser.

###  Help

To see advanced options, run
```bash
docker run -t evolbioinfo/pastml -h
```
