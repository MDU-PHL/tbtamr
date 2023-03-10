import pytest,pathlib, logging,collections,json
from unittest.mock import patch, PropertyMock

from tbtamr.AmrSetup import AmrSetup
from tbtamr.RunProfiler import RunProfiler
from tbtamr.Collate import Inferrence

test_folder = pathlib.Path(__file__).parent.parent / 'tests'
tbp_res = json.load(open(f"{test_folder}/results/tb-profiler_report.json"))

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
        print(args)
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
        # print(args)
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj._cwd = pathlib.Path(__file__).parent.parent

        assert amr_obj._check_output_file( seq_id = 'tests', step = 'profile') == f'{test_folder}/results/tests.results.json'

Header_MDU = header =  [
                "Seq_ID",
                "Species",
                "Phylogenetic lineage",
                'Predicted drug resistance',
                "Rifampicin",
                "Rifampicin - Interpretation",
                "Rifampicin - Confidence",
                "Isoniazid",
                "Isoniazid - Interpretation",
                "Isoniazid - Confidence",
                "Pyrazinamide",
                "Pyrazinamide - Interpretation",
                "Pyrazinamide - Confidence",
                "Ethambutol",
                "Ethambutol - Interpretation",
                "Ethambutol - Confidence",
                "Moxifloxacin",
                "Moxifloxacin - Interpretation",
                "Moxifloxacin - Confidence",
                "Amikacin",
                "Amikacin - Interpretation",
                "Amikacin - Confidence",
                "Cycloserine",
                "Cycloserine - Interpretation",
                "Cycloserine - Confidence",
                "Ethionamide",
                "Ethionamide - Interpretation",
                "Ethionamide - Confidence",
                "Para-aminosalicylic acid",
                "Para-aminosalicylic acid - Interpretation",
                "Para-aminosalicylic acid - Confidence",
                "Kanamycin",
                "Kanamycin - Interpretation",
                "Kanamycin - Confidence",
                "Streptomycin",
                "Streptomycin - Interpretation",
                "Streptomycin - Confidence",
                "Capreomycin",
                "Capreomycin - Interpretation",
                "Capreomycin - Confidence",
                "Clofazimine",
                "Clofazimine - Interpretation",
                "Clofazimine - Confidence",
                "Delamanid",
                "Delamanid - Interpretation",
                "Delamanid - Confidence",
                "Bedaquiline",
                "Bedaquiline - Interpretation",
                "Bedaquiline - Confidence",
                "Linezolid",
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

        _dict = [{"Seq_ID": "test_1", 
                  "Rifampicin": "No mechanism identified", 
                  "Isoniazid": "No mechanism identified", 
                  "Pyrazinamide": "No mechanism identified", 
                  "Ethambutol": "No mechanism identified", 
                  "Moxifloxacin": "No mechanism identified", 
                  "Amikacin": "No mechanism identified", 
                  "Cycloserine": "No mechanism identified", 
                  "Ethionamide": "No mechanism identified", 
                  "Para-aminosalicylic acid": "No mechanism identified", 
                  "Clofazimine": "No mechanism identified", 
                  "Delamanid": "No mechanism identified", 
                  "Bedaquiline": "No mechanism identified", 
                  "Linezolid": "No mechanism identified", 
                  "Kanamycin": "No mechanism identified", 
                  "Streptomycin": "No mechanism identified", 
                  "Capreomycin": "No mechanism identified", 
                  "Predicted drug resistance": "No drug resistance predicted", 
                  "Species": "Not likely M. tuberculosis", 
                  "Phylogenetic lineage": "Not typable", 
                  "Database version": "No database version available", 
                  "Median genome coverage": 98, "Percentage reads mapped": 98.7}]
        
        assert amr_obj._open_json(path=f'{test_folder}/results/data/test_1/tbtamr.json') == _dict

def test_get_prop_mtb_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='sample_id',res = tbp_res,val = 'median_coverage') == 98

def test_get_prop_mtb_fail():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='sample_id',res = tbp_res,val = 'median_coverage') != 97

        
def test_get_cov_success():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 40, 80) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._get_qc_feature(seq_id='sample_id',res = tbp_res,val = 'pct_reads_mapped') == 98.7

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

        assert amr_obj._check_quality(cov = 35, perc = 90) == 'Fail QC'


def test_get_qual_fail_low_prop_mtb():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 50, perc = 79) == 'Fail QC'

def test_get_qual_fail_both():

    with patch.object(AmrSetup, "__init__", lambda x: None):
        to_input = Input(['test'], False, 80, 40) 
        amr_obj = Inferrence(to_input)
        amr_obj.logger = logging.getLogger(__name__)

        assert amr_obj._check_quality(cov = 35, perc = 79) == 'Fail QC'

#         assert amr_obj._get_qc_feature(seq_id='sample_id',res = tbp_res,val = 'pct_reads_mapped') != 100
