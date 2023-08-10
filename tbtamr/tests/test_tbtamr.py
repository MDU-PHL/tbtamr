import pytest,pathlib, logging,collections,json
from unittest import TestCase
from unittest.mock import patch, PropertyMock
import string
from tbtamr.AmrSetup import AmrSetup
from tbtamr.RunProfiler import RunProfiler
from tbtamr.Collate import Inferrence,Parse
from tbtamr.TbTamr import Tbtamr
from tbtamr.version import db_version

test_folder = pathlib.Path(__file__).parent.parent / 'tests'
tbp_res = json.load(open(f"{test_folder}/results/data/A/tb-profiler_report.json"))
data_folder = test_folder / 'results' /'data'
data_folder_single  = data_folder / 'A'
data_folder_batch = data_folder / 'iso.list'
isos_list = [i for i in string.ascii_uppercase if i not in ['J','Z']]
# isos_list = ["D"]

with open(f"{test_folder / 'tbtamr_truth.json'}", "r") as j:
    batch_truth = json.load(j)


def test_file_present():
    """
    assert true when the input file is true
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "dummy_R1.fq.gz"
    
        assert amr_obj._file_present(p)

def test_file_not_present():
    """
    assert true when the input file is true
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "dummy_R5.fq.gz"
    
        assert not amr_obj._file_present(p)


def test_resorces_maxjob_1():
    """
    assert true when the maxjob is returned as the number of tbprofiler because jobs is not provided
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.one = 16
        amr_obj.five = 24
        amr_obj.fifteen = 32
        amr_obj.total_cores = 64
        val = 0

        assert amr_obj._set_threads(jobs=val) == amr_obj.fifteen/2

def test_resorces_maxjob_2():
    """
    assert true when the maxjob is returned as the number of tbprofiler because jobs is greater than available resources
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.one = 16
        amr_obj.five = 24
        amr_obj.fifteen = 32
        amr_obj.total_cores = 64
        val = 8

        assert amr_obj._set_threads(jobs=val) == int(val/2)

def test_resorces_job():
    """
    assert true when the job is returned as the number of tbprofiler because jobs is less than available resources
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.one = 16
        amr_obj.five = 24
        amr_obj.fifteen = 32
        amr_obj.total_cores = 64
        val = 3

        assert amr_obj._set_threads(jobs=val) == int(val/2)

def test_run_cmd_success():
    """
    assert true when the run cmd has exitcode == 0
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        cmd = "echo Hello World"

        assert amr_obj._run_cmd(cmd = cmd)

def test_run_cmd_failure():
    """
    assert true when the run cmd has exitcode != 0
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        cmd = "ks"

        with pytest.raises(SystemExit):
            amr_obj._run_cmd(cmd = cmd)

def test_check_prefix_success():
    """
    assert true when the prefix is provided
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.prefix = 'some_prefix'
        assert amr_obj._check_prefix()

def test_check_prefix_success():
    """
    raise system exit the prefix is not
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.prefix = ''
        with pytest.raises(SystemExit):
            amr_obj._check_prefix()

def test_reads_exist_success():
    """
    assert true when reads in input file are present
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        read = test_folder / "dummy_R1.fq.gz"

        
        assert amr_obj._check_reads(read=read)

def test_reads_exist_fail():
    """
    assert true when reads in input file are present
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        amr_obj = AmrSetup()
        amr_obj.logger = logging.getLogger(__name__)
        read = test_folder / "dummy_R5.fq.gz"

        with pytest.raises(SystemExit):
            amr_obj._check_reads(read=read)
     


DATA = collections.namedtuple('Data', ['input_data', 'jobs', 'db', 'keep', 'keep_bam','exclude_not_reportable', 'min_depth', 'min_cov','prop_mtb'])


def test_generate_cmd_batch_success():
    """
    assert True when non-empty string is given
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        args = DATA(f"{test_folder / 'isolates.tab'}", 3, '--db tbdb',False,False,False,20, 40,80)
        amr_obj = RunProfiler(args)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._batch_cmd(input_data=args.input_data) == f"parallel --colsep '\\t' -j {args.jobs} 'tb-profiler profile --read1 {{2}} --read2 {{3}} {args.db} --prefix {{1}} --dir {{1}} --min_depth 20 --no_trim --call_whole_genome --threads 1 >> {{1}}/tbprofiler.log 2>&1' :::: {args.input_data}"


# testing collation
Input = collections.namedtuple('Input', ['isolates',  'exclude_not_reportable', 'prop_mtb', 'min_cov'])
                
def test_check_output_file_profile_success():
    """
    assert True when the tb-profiler output is present
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj._cwd = pathlib.Path(__file__).parent.parent

        assert amr_obj._check_output_file( seq_id = 'tests', step = 'profile') == f'{test_folder}/results/tests.results.json'

Header_MDU = header =  [
                "Seq_ID",
                "Species",
                "Phylogenetic lineage",
                'Predicted drug resistance',
                "Rifampicin - ResMech",
                "Rifampicin - Interpretation",
                "Rifampicin - Confidence",
                "Isoniazid - ResMech",
                "Isoniazid - Interpretation",
                "Isoniazid - Confidence",
                "Pyrazinamide - ResMech",
                "Pyrazinamide - Interpretation",
                "Pyrazinamide - Confidence",
                "Ethambutol - ResMech",
                "Ethambutol - Interpretation",
                "Ethambutol - Confidence",
                "Moxifloxacin - ResMech",
                "Moxifloxacin - Interpretation",
                "Moxifloxacin - Confidence",
                "Amikacin - ResMech",
                "Amikacin - Interpretation",
                "Amikacin - Confidence",
                "Cycloserine - ResMech",
                "Cycloserine - Interpretation",
                "Cycloserine - Confidence",
                "Ethionamide - ResMech",
                "Ethionamide - Interpretation",
                "Ethionamide - Confidence",
                "Para-aminosalicylic acid - ResMech",
                "Para-aminosalicylic acid - Interpretation",
                "Para-aminosalicylic acid - Confidence",
                "Kanamycin - ResMech",
                "Kanamycin - Interpretation",
                "Kanamycin - Confidence",
                "Streptomycin - ResMech",
                "Streptomycin - Interpretation",
                "Streptomycin - Confidence",
                "Capreomycin - ResMech",
                "Capreomycin - Interpretation",
                "Capreomycin - Confidence",
                "Clofazimine - ResMech",
                "Clofazimine - Interpretation",
                "Clofazimine - Confidence",
                "Delamanid - ResMech",
                "Delamanid - Interpretation",
                "Delamanid - Confidence",
                "Bedaquiline - ResMech",
                "Bedaquiline - Interpretation",
                "Bedaquiline - Confidence",
                "Linezolid - ResMech",
                "Linezolid - Interpretation",
                "Linezolid - Confidence",
                "Quality",
                "Database version",
            ]

Header_General = header =  [
                "Seq_ID",
                "Species",
                "Phylogenetic lineage",
                'Predicted drug resistance',
                "Rifampicin",
                "Isoniazid",
                "Pyrazinamide",
                "Ethambutol",
                "Moxifloxacin",
                "Amikacin",
                "Cycloserine",
                "Ethionamide",
                "Kanamycin",
                "Streptomycin",
                "Capreomycin",
                "Para-aminosalicylic acid",
                "Clofazimine",
                "Delamanid",
                "Bedaquiline",
                "Linezolid",
                "Median genome coverage",
                "Percentage reads mapped",
                "Quality",
                "Database version"
            ]

def test_check_headers_mdu_success():
    """
    assert True when headers are correct
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_header(_type = 'mdu') == Header_MDU

def test_check_headers_mdu_fail():
    """
    assert True when headers are correct
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_header(_type = 'general') != Header_MDU

def test_check_headers_mdu_fail():
    """
    assert True when headers are correct
    """
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_header(_type = 'general') == Header_General
    
def test_open_json_fail():
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        amr_obj._cwd = pathlib.Path(__file__).parent.parent

        with pytest.raises(SystemExit):
            amr_obj._open_json(path=f'{test_folder}/results/tests.results.json')


def test_open_json_success():
    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        amr_obj._cwd = pathlib.Path(__file__).parent.parent

        _dict = {"A": 
                 {"rifampicin": "-", 
                  "isoniazid": "-", 
                  "pyrazinamide": "-", 
                  "ethambutol": "-", 
                  "streptomycin": "-", 
                  "fluoroquinolones": "-", 
                  "moxifloxacin": "-", 
                  "ofloxacin": "-", 
                  "levofloxacin": "-", 
                  "ciprofloxacin": "-", 
                  "aminoglycosides": "-", 
                  "amikacin": "-", 
                  "kanamycin": "-", 
                  "capreomycin": "-", 
                  "ethionamide": "-", 
                  "para-aminosalicylic_acid": "-", 
                  "cycloserine": "-", 
                  "linezolid": "-", 
                  "bedaqualine": "-", 
                  "bedaquiline": "-", 
                  "clofazamine": "-", 
                  "clofazimine": "-", 
                  "delamanid": 
                  "fbiA_p.Thr302Met", 
                  "main_lin": 
                  "lineage4", 
                  "sublin": "lineage4.1.1.1", 
                  "drtype": "Other", 
                  "pct_reads_mapped": 98.21, 
                  "num_reads_mapped": 4262601, 
                  "median_coverage": 145, 
                  "num_dr_variants": 1, 
                  "num_other_variants": 18}
                  }
        assert amr_obj._open_json(path=f'{test_folder}/results/data/A/tb-profiler_report.json') == _dict

def test_get_prop_mtb_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='A',res = tbp_res,val = 'median_coverage') == 145

def test_get_prop_mtb_fail():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='A',res = tbp_res,val = 'median_coverage') != 97

        
def test_get_cov_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='A',res = tbp_res,val = 'pct_reads_mapped') == 98.21

def test_get_qual_pass():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 50, perc = 90) == 'Pass QC'

        

def test_get_qual_fail_low_cov():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 35, perc = 90) == 'Failed: < 40x aligned coverage to reference genome'


def test_get_qual_fail_low_prop_mtb():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 50, perc = 79) == 'Failed: < 80 % _M. tuberculosis_ reads in sample'

def test_get_qual_fail_both():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 35, perc = 79) == 'Failed: < 80 % _M. tuberculosis_ reads in sample'

def test_check_mutant_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "ethambutol"
        mut = "embB_p.Leu74Arg"

        assert amr_obj._check_mut(drug = drug, mut=mut) == "embB_p.Leu74Arg"


def test_check_mutant_success_large_deletion():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "isoniazid"
        mut = "katG_del"

        assert amr_obj._check_mut(drug = drug, mut=mut) == "katG_large_deletion"


def test_check_mutant_fail():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "rifampicin"
        mut = "rpoB_del"

        assert amr_obj._check_mut(drug = drug, mut=mut) == "No mechanism identified*"


def test_check_confidence_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "ethambutol"
        mut = "embB_p.Leu74Arg"

        assert amr_obj._get_confidence(drug = drug, mut=mut,qual = "Pass qc") == "Moderate"

def test_check_confidence_fail_qc():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "ethambutol"
        mut = "embB_p.Leu74Arg"

        assert amr_obj._get_confidence(drug = drug, mut=mut,qual = "Fail qc") == ""


def test_check_interpretation_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "ethambutol"
        mut = "embB_p.Leu74Arg"

        assert amr_obj._get_interpret(drug = drug, mut=mut,qual = "Pass qc") == "Resistant"



def test_check_interpretation_fail_qc():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        drug = "ethambutol"
        mut = "embB_p.Leu74Arg"

        assert amr_obj._get_interpret(drug = drug, mut=mut,qual = "Fail qc") == ""


def test_check_output_file_collate_success():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        seq_id = "A"
        step = "collate"


        assert amr_obj._check_output_file(seq_id=seq_id, step = step) == f"{data_folder_single}/tb-profiler_report.json"



def test_check_output_file_collate_fail():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        seq_id = "Z"
        step = "collate"

        assert amr_obj._check_output_file(seq_id=seq_id, step = step) == False



def test_check_output_file_profile_success():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        seq_id = "A"
        step = "profile"


        assert amr_obj._check_output_file(seq_id=seq_id, step = step) == f"{data_folder_single}/results/A.results.json"



def test_check_output_file_collate_fail():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        seq_id = "Z"
        step = "profile"

        assert amr_obj._check_output_file(seq_id=seq_id, step = step) == False

def test_check_output_profile():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        isolates = {"A":{}}
        step = "profile"

        assert amr_obj._check_output(isolates=isolates, step = step) == {"A":{'profile':f"{data_folder_single}/results/A.results.json"}}


def test_check_output_collate():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        
        amr_obj = Tbtamr()
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        isolates = {"A":{}}
        step = "collate"

        assert amr_obj._check_output(isolates=isolates, step = step) == {"A":{'collate':f"{data_folder_single}/tb-profiler_report.json"}}



def test_fail_isolate_dict():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input('sample_id', False, 40, 80) 
        amr_obj = Parse(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        with pytest.raises(SystemExit):
            amr_obj._fail_isolate_dict(seq_id='sample_id',step = 'collate')




def test_get_isolate_dict_success():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        to_input = Input('sample_id', False, 40, 80) 
        amr_obj = Parse(to_input)
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        isos = ["A"]
        # step = "collate"

        assert amr_obj._get_isolate_dict(isos=isos) == {"A":{'collate':f"{data_folder_single}/tb-profiler_report.json", 'profile':f"{data_folder_single}/results/A.results.json"}}


def test_get_isolate_dict_fail():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        to_input = Input('sample_id', False, 40, 80) 
        amr_obj = Parse(to_input)
        amr_obj._cwd = data_folder
        amr_obj.logger = logging.getLogger(__name__)
        # amr_obj.db = pathlib.Path(__file__).parent.parent /"db" /
        isos = ["Z"]
        # step = "collate"
        with pytest.raises(SystemExit):
            assert amr_obj._get_isolate_dict(isos=isos)

def test_interpretations():

    with patch.object(Tbtamr, "__init__", lambda x: None):
        for i in isos_list:
            print(i)
            to_input = Input(i, False, 80, 40) 
            truth = batch_truth[i]
            amr_obj = Parse(to_input)
            amr_obj._cwd = data_folder
            amr_obj.logger = logging.getLogger(__name__)
            isos = [i]
            isolates = amr_obj._get_isolate_dict(isos=isos)
            to_infer = Input(isolates,False, 80,40)
            infer_obj = Inferrence(to_infer)
            infer_obj._cwd = data_folder
            db_path = f"{db_version}.json"
            infer_obj.db_path = f"{pathlib.Path(__file__).parent.parent /'db' / db_path}"
            _dict = infer_obj._infer_single_seq(isolate=i)
            res = infer_obj._wrangle_json(res = _dict)
            assert res == truth