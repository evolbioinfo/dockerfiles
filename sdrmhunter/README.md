# SDRM Hunter

This docker file wraps [sierrapy](https://github.com/hivdb/sierra-client/tree/master/python) python3 module
to detect SDRMs in a given HIV fasta sequence alignment.

### How to run SDRM Hunter
As an input, one needs to provide an HIV sequence alignment in fasta format (to be specified with *--fasta* option).

Basic usage:
```bash
docker run -v <path_to_the_folder_containing_the_alignment>:/data:rw -t evolbioinfo/sdrmhunter --fasta /data/<alignment_file> --output /data/<sdrm_annotation_file>
```

###  Help

To see the options, run
```bash
docker run -t evolbioinfo/sdrmhunter -h
```