import pytest,pathlib, logging,collections,json,string
from unittest import TestCase
from unittest.mock import patch, PropertyMock

from tbtamr.Predict import PredictAmr

test_folder = pathlib.Path(__file__).parent.parent / 'tests'

def test_check_file():
    """
    assert true when the input file is true
    """
    with patch.object(PredictAmr, "__init__", lambda x: None):
        amr_obj = PredictAmr()
        amr_obj.logger = logging.getLogger(__name__)
        p = test_folder / "empty_file.txt"
    
        assert amr_obj._file_present(p)

