import argparse, sys, pathlib, tempfile
from distutils.command.install_egg_info import to_filename
from Parse import Vcf
from Predict import PredictAmr
from Call import generatevcf

# from tbtamr.RunProfiler import RunProfiler
# from tbtamr.Collate import Inferrence, Parse, Mdu
# from tbtamr.TbTamr_Utils import check,install
from version import __version__, db_version

"""
tbtamr is designed to implement TB-profiler and parse the results compatible for MDU use. It may also be used for other purposes where the format of output is compatible

"""

def run_predict(args):
    if args.vcf == '':
        vcf = generatevcf(read1 = args.read1,
                        read2 = args.read2,
                        threads = args.threads,
                        ram = args.ram,
                        seq_id = args.seq_id,
                        keep = args.keep,
                        mindepth = args.min_depth,
                        minfrac = args.min_frac,
                        force = args.force,
                        mtb = True,
                        tmp = args.tmp
                        )
    else:
        vcf = args.vcf
    Prs = Vcf(vcf = vcf,
              catalog= args.catalog,
              catalog_config = args.catalog_config,
              seq_id=args.seq_id,
              )
    variants = Prs.get_variant_data()

    P = PredictAmr(variants = variants,
                 catalog = args.catalog,
                 config = args.catalog_config,
                 rules = args.rules,
                 seq_id = args.seq_id,
                 vcf = args.vcf,
                 ref = args.reference_file,
                 barcode = args.barcode
                )
    P.run_prediction()

def search_catalog(args):
    pass

def set_parsers():
    parser = argparse.ArgumentParser(
        description="Genomic AMR prediction pipeline for Mtb for clinical and public health reporting", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    
    subparsers = parser.add_subparsers(help="Task to perform")
    
    parser_sub_predict = subparsers.add_parser('predict', help='Predict AMR results', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # parser_sub_predict.add_argument()
    
    parser_sub_predict.add_argument(
        '--vcf',
        help="VCF file generated using the H37rV v3 reference genome",
        default=""
    )
    parser_sub_predict.add_argument(
        '--read1',
        '-1',
        help="Path to read 1 - not required if vcf is supplied",
        default=""
    )
    parser_sub_predict.add_argument(
        '--read2',
        '-2',
        help="Path to read 2 - not required if vcf is supplied",
        default=""
    )
    parser_sub_predict.add_argument(
        "--catalog",
        "-c",
        help="csv variant catalog",
        # required=True,
        default=f"{pathlib.Path(__file__).parent / 'db'/ 'who_v2_catalog.csv'}"
    )
    parser_sub_predict.add_argument(
        '--catalog_config',
        '-cfg',
        # required=True,
        help = "json file indicating the relevant column settings for interpretation of the catalog file.",
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'db_config.json'}"
    )
    parser_sub_predict.add_argument(
        '--rules',
        '-r',
        # required = True,
        help= f"csv file with rules for predicting resistance profiles from genomic data.",
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'rules.csv'}"
    )
    parser_sub_predict.add_argument(
        '--barcode',
        '-b',
        # required = True,
        help= f"Barcode to use for lineage calling and speciation.",
        default= f"{pathlib.Path(__file__).parent / 'db'/ 'tbtamr.barcode.bed'}"
    )
    parser_sub_predict.add_argument(
        '--reference_file',
        '-ref',
        # required = True,
        help= f"Reference file to use for calling lineage.",
        default= f"{pathlib.Path(__file__).parent / 'db'/ 'tbtamr.fasta'}"
    )
    parser_sub_predict.add_argument(
        '--seq_id',
        '-s',
        # required=True,
        help= "Sequence name.",
        default=""
    )
    parser_sub_predict.add_argument(
        '--min_depth',
        '-md',
        help= f"Minimum depth to call a variant - only required if reads supplied as input",
        default= 20
    )
    parser_sub_predict.add_argument(
        '--min_frac',
        '-mf',
        help= f"Minimum proportion to call a variant (0-1) - only required if reads supplied as input",
        default= 0.1
    )
    parser_sub_predict.add_argument(
        '--threads',
        '-t',
        help = "Threads to use for generation of vcf file - only required if reads supplied as input",
        default = 8
    )
    parser_sub_predict.add_argument(
        '--ram',
        help = "Max ram to use - only required if reads supplied as input",
        default = 8
    )
    parser_sub_predict.add_argument(
        '--keep',
        '-k',
        help = "Keep accessory files generated from alignment for further use.",
        action = "store_true"
    )
    parser_sub_predict.add_argument(
        '--force',
        '-f',
        help = "Force replace an existing folder.",
        action = "store_true"
    )
    parser_sub_predict.add_argument(
        '--cascade',
        action='store_true',
        help = 'If you would like to apply cascade reporting structure.'
    )
    parser_sub_predict.add_argument(
        '--report_rules',
        help = 'Path to file describing reporting rules.',
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'report_rules.csv'}"
    )
    parser_sub_predict.add_argument(
        '--tmp',
        help = "temp directory to use",
        default = f"{pathlib.Path(tempfile.gettempdir())}"
    )
    # parser_sub_search = subparsers.add_parser('search', help='Search the provided catalog for variant information', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser_sub_search.add_argument(
    #     "--catalog",
    #     "-c",
    #     help="csv variant catalog",
    #     # required=True,
    #     default=f"{pathlib.Path(__file__).parent / 'db'/ 'who_v2_catalog.csv'}"
    # )
    # parser_sub_search.add_argument(
    #     '--catalog_config',
    #     '-cfg',
    #     # required=True,
    #     help = "json file indicating the relevant column settings for interpretation of the catalog file.",
    #     default= f"{pathlib.Path(__file__).parent / 'configs'/ 'db_config.json'}"
    # )
    # parser_sub_search.add_argument(
    #     '--query',
    #     '-q',
    #     required=True,
    #     nargs='+',
    #     help="The term and column to search. Example rifampicin drug - this will search for rifampicin in the drug column"
    # )
    # parser_sub_search.set_defaults(func = search_catalog)
    parser_sub_predict.set_defaults(func=run_predict)
    args = parser.parse_args(args=None if sys.argv[1:]  else ['--help'])
    return args

 
def main():
    """
    run pipeline
    """

    args = set_parsers()
    args.func(args)
    

if __name__ == "__main__":
    main()
