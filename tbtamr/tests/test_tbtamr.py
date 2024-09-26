import pytest,pathlib, logging,collections,json,string,pandas
from unittest import TestCase
from unittest.mock import patch, PropertyMock

from tbtamr.Predict import PredictAmr

test_folder = pathlib.Path(__file__).parent.parent / 'tests'
rules_path = pathlib.Path(__file__).parent.parent / 'configs' 

classification_rules = pandas.read_csv(f"{rules_path / 'classification_criteria.csv'}")
interpretation_rules = pandas.read_csv(f"{test_folder / 'interpretation_criteria.csv'}")

def test_check_file():
    """
    assert true when the input file is true
    """
    with patch.object(PredictAmr, "__init__", lambda x: None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "empty_file.txt"
    
        assert amr_obj.check_file(p)

def test_check_file_fail():
    """
    assert true when the input file is not present - test proper failure
    """
    with patch.object(PredictAmr, "__init__", lambda x: None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "no_empty_file.txt"
    
        with pytest.raises(SystemExit):
            amr_obj.check_file(p)
    
def test_collect_resistance_mechs():
    """
    assert true when dataframe is appropriately filtered
    """
    with patch.object(PredictAmr, "__init__", lambda x: None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.config = {"variant_col":"variant"}
        variants = [{
            "variant":"y"
        }]
        catalog = pandas.DataFrame({"variant":"x","drug":"a"}, index = [0])

        assert amr_obj.collect_resistance_mechs(catalog=catalog, variants=variants).empty

def test_get_highest_conf_empty():
    """
    assert true when no conf is sent - return empty string
    """
    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.config = {
            "confidence_levels": {
            "1) Assoc w R" : 0,
            "2) Assoc w R - Interim": 1
            },
            "confidence_key":{
            "1) Assoc w R" : "High",
            "2) Assoc w R - Interim": "Moderate",
            "3) Uncertain significance": "Uncertain"
            }
        }
        conf = []

        assert amr_obj.get_highest_conf(conf=conf) == ""
    
def test_get_highest_conf_high():
    """
    assert true when no conf is sent - return empty string
    """
    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.config = {
            "confidence_levels": {
            "1) Assoc w R" : 0,
            "2) Assoc w R - Interim": 1
            },
            "confidence_key":{
            "1) Assoc w R" : "High",
            "2) Assoc w R - Interim": "Moderate",
            "3) Uncertain significance": "Uncertain"
            }
        }
        conf = ["1) Assoc w R","2) Assoc w R - Interim"]

        assert amr_obj.get_highest_conf(conf=conf) == "High"

    
def test_get_highest_conf_mod():
    """
    assert true when no conf is sent - return empty string
    """
    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.config = {
            "confidence_levels": {
            "1) Assoc w R" : 0,
            "2) Assoc w R - Interim": 1
            },
            "confidence_key":{
            "1) Assoc w R" : "High",
            "2) Assoc w R - Interim": "Moderate",
            "3) Uncertain significance": "Uncertain"
            }
        }
        conf = ["2) Assoc w R - Interim"]

        assert amr_obj.get_highest_conf(conf=conf) == "Moderate"

def test_get_rules_for_dr():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        interp_truth = interpretation_rules[interpretation_rules['drug'] == 'amikacin']
        interp_truth = interp_truth.fillna('')
        dr = 'amikacin'
        assert amr_obj.get_rules_for_dr(dr = dr, rules = interpretation_rules).equals(interp_truth)
    

def test_get_rules_for_dr_fail():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        interp_truth = interpretation_rules[interpretation_rules['drug'] == 'amikacin']
        interp_truth = interp_truth.fillna('')
        dr = 'rifampicin'
        assert not amr_obj.get_rules_for_dr(dr = dr, rules = interpretation_rules).equals(interp_truth)

def test_check_shape_no_shape():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        rule = ""

        assert not amr_obj.check_shape(rule)


def test_check_shape_has_shape():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        rule = "==1"

        assert amr_obj.check_shape(rule)

def test_extract_mutations_multiple():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        amr_obj.config = {"confidence_levels": {
                            "1) Assoc w R" : 0,
                            "2) Assoc w R - Interim": 1
            },
            "confidence_key":{
                    "1) Assoc w R" : "High",
                    "2) Assoc w R - Interim": "Moderate",
                    "3) Uncertain significance": "Uncertain"
                    }
            }
        result = {"rifampicin - mechanisms": "var1 (High);var2 (Moderate)"}

        assert amr_obj.extract_mutations(result=result, dr = 'rifampicin') == ["var1","var2"]


def test_extract_mutations_none():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        result = {"rifampicin - mechanisms": 'No reportable mechanims'}

        assert amr_obj.extract_mutations(result=result, dr = 'rifampicin') == []

def test_update_result():

    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        result = {"rifampicin - interpretation": "Resistant"}
        row = 0,interpretation_rules.iloc[6]
        
        assert amr_obj.update_result(dr = "rifampicin", result=result, rule = row) == {"rifampicin - interpretation":"Low-level resistant", "rifampicin - override": "if there is a single variant that equals rpoB_p.Asp435Tyr then the interpretation is low-level resistant - this will override the default rule - custom"}


def test_construct_rule():
    with patch.object(PredictAmr, "__init__", lambda x:None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        row = 0,interpretation_rules.iloc[6]

        assert amr_obj.construct_rule(row = row) == "`variant` in ['rpoB_p.Asp435Tyr']"
