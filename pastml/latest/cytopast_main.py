from cytopast.pastml_analyser import pastml_pipeline

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser(description="Visualisation of annotated phylogenetic trees (as html maps).")

    annotation_group = parser.add_argument_group('annotation-related arguments')
    annotation_group.add_argument('--data', required=True, type=str,
                                  help="the annotation file in tab/csv format with the first row "
                                       "containing the column names.")
    annotation_group.add_argument('--data_sep', required=False, type=str, default='\t',
                                  help="the column separator for the data table. "
                                       "By default is set to tab, i.e. for tab file. " \
                                       "Set it to ',' if your file is csv.")
    annotation_group.add_argument('--id_index', required=False, type=int, default=0,
                                  help="the index of the column in the data table that contains the tree tip names, "
                                       "indices start from zero (by default is set to 0).")
    annotation_group.add_argument('--columns', nargs='*',
                                  help="names of the data table columns that contain states "
                                       "to be analysed with PASTML. "
                                       "If neither columns nor copy_columns are specified, "
                                       "then all columns will be considered for PASTMl analysis.",
                                  type=str)
    annotation_group.add_argument('--copy_columns', nargs='*',
                                  help="names of the data table columns that contain states to be copied as-is, "
                                       "without applying PASTML (the missing states will stay unresolved).",
                                  type=str)

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('--tree', help="the input tree in newick format.", type=str, required=True)

    pastml_group = parser.add_argument_group('ancestral-state inference-related arguments')
    pastml_group.add_argument('--model', required=False, default='JC', type=str,
                              help="the evolutionary model to be used by PASTML (can be JC or F81).")
    pastml_group.add_argument('--work_dir', required=False, default=None, type=str,
                              help="the working dir for PASTML to put intermediate files into "
                                   "(if not specified a temporary dir will be created).")
    pastml_group.add_argument('--cache', action='store_true',
                              help="if set, the results of previous PASTML runs on this data will be reused "
                                   "when possible")

    vis_group = parser.add_argument_group('visualisation-related arguments')
    vis_group.add_argument('--name_column', type=str, default=None,
                           help="name of the data table column to be used for node names "
                                "in the compressed map visualisation"
                                "(must be one of those specified in columns or copy_columns if they are specified)."
                                "If the data table contains only one column it will be used by default.")
    vis_group.add_argument('--all', action='store_true', help="Keep all the nodes in the compressed map visualisation, "
                                                              "even the minor ones.")

    out_group = parser.add_argument_group('output-related arguments')
    out_group.add_argument('--out_data', required=False, type=str,
                           help="the output annotation file with the states inferred by PASTML.")
    out_group.add_argument('--html_compressed', required=False, default=None, type=str,
                           help="the output summary map visualisation file (html).")
    out_group.add_argument('--html', required=False, default=None, type=str,
                           help="the output tree visualisation file (html).")

    parser.add_argument('--verbose', action='store_true', help="print information on the progress of the analysis")
    params = parser.parse_args()

    pastml_pipeline(**vars(params))

