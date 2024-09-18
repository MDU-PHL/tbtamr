import argparse, sys, pathlib, tempfile, os
from distutils.command.install_egg_info import to_filename
from Parse import Vcf
from Predict import PredictAmr
from Utils import check_annotate, check_mutamr, check_lineage
from Annotate import annotate
from Search import search
from version import __version__, db_version

"""
tbtAMR implements user defined (or default) criteria for the inference of phenotypic AMR in M. tuberculosis.

"""
def run_predict(args):

    Prs = Vcf(vcf = args.vcf,
              catalog= args.catalog,
              catalog_config = args.catalog_config,
              seq_id=args.seq_id,
              )
    variants = Prs.get_variant_data()
    call_lineage = False
    if args.call_lineage and check_lineage():
        call_lineage = True
    P = PredictAmr(variants = variants,
                 catalog = args.catalog,
                 config = args.catalog_config,
                 interpretation_rules = args.interpretation_rules,
                 classification_rules = args.classification_rules,
                 seq_id = args.seq_id,
                 vcf = args.vcf,
                 ref = args.reference_file,
                 barcode = args.barcode,
                 cascade = args.cascade,
                 call_lineage = call_lineage
                )
    P.run_prediction()


def run_fq2vcf(args):
    
    from Call import generatevcf

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

def run_annotate(args):

    vcf_file = annotate(vcf_file = args.vcf,
                         seq_id= args.seq_id)

def run_full(args):
    from Call import generatevcf
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
    call_lineage = False
    if args.call_lineage and check_lineage():
        call_lineage = True
    P = PredictAmr(variants = variants,
                 catalog = args.catalog,
                 config = args.catalog_config,
                 interpretation_rules = args.interpretation_rules,
                 classification_rules = args.classification_rules,
                 seq_id = args.seq_id,
                 vcf = vcf,
                 ref = args.reference_file,
                 barcode = args.barcode,
                 cascade = args.cascade,
                 call_lineage = call_lineage
                )
    P.run_prediction()

def search_catalog(args):
    
    search(config=args.catalog_config,
            catalogue=args.catalog,
            query=args.query)



def set_parsers():
    parser = argparse.ArgumentParser(
        description="Genomic AMR prediction pipeline for Mtb for clinical and public health reporting", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    
    subparsers = parser.add_subparsers(help="Task to perform")
    parser_sub_search = subparsers.add_parser('search', help='Search the provided catalog for variant information', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub_search.add_argument(
        "--catalog",
        "-c",
        help="csv variant catalog",
        # required=True,
        default=f"{pathlib.Path(__file__).parent / 'db'/ 'who_v2_catalog.csv'}"
    )
    parser_sub_search.add_argument(
        '--catalog_config',
        '-cfg',
        # required=True,
        help = "json file indicating the relevant column settings for interpretation of the catalog file.",
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'db_config.json'}"
    )
    parser_sub_search.add_argument(
        '--query',
        '-q',
        required=True,
        nargs='+',
        help="The variant to search for. Multiple allowed separated by a space."
    )
    parser_sub_predict = subparsers.add_parser('predict', help='Predict AMR results from vcf', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            # parser_sub_predict.add_argument()
    parser_sub_predict.add_argument(
        '--seq_id',
        '-s',
        # required=True,
        help= "Sequence name.",
        default="tbtamr"
    )
    parser_sub_predict.add_argument(
        '--vcf',
        help="VCF file generated using the H37rV v3 reference genome",
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
        '--interpretation_rules',
        '-r',
        # required = True,
        help= f"csv file with rules for predicting resistance profiles from genomic data.",
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'interpretation_rules.csv'}"
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
        '--classification_rules',
        '-cr',
        # required = True,
        help= f"csv file with rules for predicting resistance profiles from genomic data.",
        default= f"{pathlib.Path(__file__).parent / 'configs'/ 'classification_rules.csv'}"
    )
    parser_sub_predict.add_argument(
        '--force',
        '-f',
        help = "Force replace an existing folder.",
        action = "store_true"
    )
    parser_sub_predict.add_argument(
        '--call_lineage',
        help = "Use pathogen profiler to call lineage",
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
    if check_annotate():
        parser_sub_annotate = subparsers.add_parser('annotate', help='Annotate a vcf file to be compatible with WHO catalogue v2 - NO report generated.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser_sub_annotate.add_argument(
                '--seq_id',
                '-s',
                # required=True,
                help= "Sequence name.",
                default="tbtamr"
            )
        parser_sub_annotate.add_argument(
            '--vcf',
            help="VCF file generated using the H37rV v3 reference genome",
            default=""
        )

        if check_mutamr():
            from Call import generatevcf
            parser_sub_fqtovcf = subparsers.add_parser('fq2vcf', help='Generate an annotated vcf file from paired-end fastq files - NO report generated.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser_sub_fqtovcf.add_argument(
                '--seq_id',
                '-s',
                # required=True,
                help= "Sequence name.",
                default="tbtamr"
            )
            parser_sub_fqtovcf.add_argument(
                '--read1',
                '-1',
                help="Path to read 1 - not required if vcf is supplied",
                default=""
            )
            parser_sub_fqtovcf.add_argument(
                '--read2',
                '-2',
                help="Path to read 2 - not required if vcf is supplied",
                default=""
            )
            parser_sub_fqtovcf.add_argument(
                '--min_depth',
                '-md',
                help= f"Minimum depth to call a variant - only required if reads supplied as input",
                default= 20
            )
            parser_sub_fqtovcf.add_argument(
                '--min_frac',
                '-mf',
                help= f"Minimum proportion to call a variant (0-1) - only required if reads supplied as input",
                default= 0.1
            )
            parser_sub_fqtovcf.add_argument(
                '--threads',
                '-t',
                help = "Threads to use for generation of vcf file - only required if reads supplied as input",
                default = 8
            )
            parser_sub_fqtovcf.add_argument(
                '--ram',
                help = "Max ram to use - only required if reads supplied as input",
                default = 8
            )
            parser_sub_fqtovcf.add_argument(
                '--keep',
                '-k',
                help = "Keep accessory files generated from alignment for further use.",
                action = "store_true"
            )
            parser_sub_fqtovcf.add_argument(
                '--force',
                '-f',
                help = "Force replace an existing folder.",
                action = "store_true"
            )
            parser_sub_fqtovcf.add_argument(
                '--tmp',
                help = "temp directory to use",
                default = f"{pathlib.Path(tempfile.gettempdir())}"
            )

            parser_sub_full = subparsers.add_parser('full', help='Predict AMR results from reads or vcf', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            # parser_sub_predict.add_argument()
            parser_sub_full.add_argument(
                '--seq_id',
                '-s',
                # required=True,
                help= "Sequence name.",
                default="tbtamr"
            )
            parser_sub_full.add_argument(
                '--vcf',
                help="VCF file generated using the H37rV v3 reference genome",
                default=""
            )
            parser_sub_full.add_argument(
                '--read1',
                '-1',
                help="Path to read 1 - not required if vcf is supplied",
                default=""
            )
            parser_sub_full.add_argument(
                '--read2',
                '-2',
                help="Path to read 2 - not required if vcf is supplied",
                default=""
            )
            parser_sub_full.add_argument(
                "--catalog",
                "-c",
                help="csv variant catalog",
                # required=True,
                default=f"{pathlib.Path(__file__).parent / 'db'/ 'who_v2_catalog.csv'}"
            )
            parser_sub_full.add_argument(
                '--catalog_config',
                '-cfg',
                # required=True,
                help = "json file indicating the relevant column settings for interpretation of the catalog file.",
                default= f"{pathlib.Path(__file__).parent / 'configs'/ 'db_config.json'}"
            )
            parser_sub_full.add_argument(
                '--interpretation_rules',
                '-r',
                # required = True,
                help= f"csv file with rules for predicting resistance profiles from genomic data.",
                default= f"{pathlib.Path(__file__).parent / 'configs'/ 'interpretation_rules.csv'}"
            )
            parser_sub_full.add_argument(
                '--barcode',
                '-b',
                # required = True,
                help= f"Barcode to use for lineage calling and speciation.",
                default= f"{pathlib.Path(__file__).parent / 'db'/ 'tbtamr.barcode.bed'}"
            )
            parser_sub_full.add_argument(
                '--reference_file',
                '-ref',
                # required = True,
                help= f"Reference file to use for calling lineage.",
                default= f"{pathlib.Path(__file__).parent / 'db'/ 'tbtamr.fasta'}"
            )
            parser_sub_full.add_argument(
                '--classification_rules',
                '-cr',
                # required = True,
                help= f"csv file with rules for predicting resistance profiles from genomic data.",
                default= f"{pathlib.Path(__file__).parent / 'configs'/ 'classification_rules.csv'}"
            )
            parser_sub_full.add_argument(
                '--call_lineage',
                help = "Use pathogen profiler to call lineage",
                action = "store_true"
            )
            parser_sub_full.add_argument(
                '--min_depth',
                '-md',
                help= f"Minimum depth to call a variant - only required if reads supplied as input",
                default= 20
            )
            parser_sub_full.add_argument(
                '--min_frac',
                '-mf',
                help= f"Minimum proportion to call a variant (0-1) - only required if reads supplied as input",
                default= 0.1
            )
            parser_sub_full.add_argument(
                '--threads',
                '-t',
                help = "Threads to use for generation of vcf file - only required if reads supplied as input",
                default = 8
            )
            parser_sub_full.add_argument(
                '--ram',
                help = "Max ram to use - only required if reads supplied as input",
                default = 8
            )
            parser_sub_full.add_argument(
                '--keep',
                '-k',
                help = "Keep accessory files generated from alignment for further use.",
                action = "store_true"
            )
            parser_sub_full.add_argument(
                '--force',
                '-f',
                help = "Force replace an existing folder.",
                action = "store_true"
            )
            parser_sub_full.add_argument(
                '--cascade',
                action='store_true',
                help = 'If you would like to apply cascade reporting structure.'
            )
            parser_sub_full.add_argument(
                '--report_rules',
                help = 'Path to file describing reporting rules.',
                default= f"{pathlib.Path(__file__).parent / 'configs'/ 'report_rules.csv'}"
            )
            parser_sub_full.add_argument(
                '--tmp',
                help = "temp directory to use",
                default = f"{pathlib.Path(tempfile.gettempdir())}"
            )
            parser_sub_full.set_defaults(func=run_full)
            parser_sub_fqtovcf.set_defaults(func=run_fq2vcf)
        parser_sub_annotate.set_defaults(func=run_annotate)
    parser_sub_search.set_defaults(func = search_catalog)
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
