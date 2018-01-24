import logging

from cytopast.pastml_analyser import pastml

if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    import argparse

    parser = argparse.ArgumentParser(description="Reconstructs ancestral states with PASTML "
                                                 "and visualizes the result as an html.")

    parser.add_argument('--tree', help="the input tree in newick format.", type=str, required=True)
    parser.add_argument('--data', required=True, type=str,
                        help="the annotation file in tab/csv format with the first row containing the column names.")
    parser.add_argument('--data_sep', required=False, type=str, default='\t',
                        help="the column separator for the data table. By default is set to tab, i.e. for tab file. "
                             "Set it to , if your file is csv.")
    parser.add_argument('--id_index', required=False, type=int, default=0,
                        help="the index of the column in the data table that contains the tree tip names, "
                             "indices start from zero (by default is set to 0).")
    parser.add_argument('--columns', nargs='*',
                        help="names of the data table columns that contain states to be analysed with PASTML,"
                             "if not specified all columns will be considered.",
                        type=str)
    parser.add_argument('--name_column', type=str, default=None,
                        help="name of the data table column to be used for node names in the visualization"
                             "(must be one of those specified in columns if columns are specified)."
                             "If the data table contains only one column it will be used by default.")
    parser.add_argument('--for_names_only', action='store_true',
                        help="If specified, and if we are to analyse multiple states (specified in columns),"
                             "and the name_column is specified,"
                             "then the name_column won't be assigned a coloured section on the nodes, "
                             "but will only be shown as node names.")
    parser.add_argument('--html_compressed', required=True, default=None, type=str,
                        help="the output summary map visualisation file (html).")
    parser.add_argument('--all', action='store_true', help="Keep all the nodes in the summary map, even the minor ones.")
    parser.add_argument('--html', required=False, default=None, type=str,
                        help="the output tree visualisation file (html). "
                             "If not specified, the full tree won't be visualized")
    parser.add_argument('--model', required=False, default='JC', type=str,
                        help="the evolutionary model to be used by PASTML: "
                             "can be JC (Jukes and Cantor, 1969) or F81 (Felsenstein 1981). By default is set to JC.")
    parser.add_argument('--work_dir', required=False, default=None, type=str,
                        help="the working dir for PASTML where it will create all intermediate files"
                             "(if not specified a temporary dir will be created).")
    params = parser.parse_args()

    pastml(**vars(params))
