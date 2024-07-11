import pytest
from omniquant.calibrate import Quantifier

def test_quantifier_initialization(quantifier_instance):
    assert isinstance(quantifier_instance, Quantifier)
    assert quantifier_instance.metabolite_name is not None
    assert quantifier_instance.cal_data is not None
    assert quantifier_instance.is_int_std is not None
    assert quantifier_instance.is_int_std_conc_known is not None
    assert quantifier_instance.is_cal_points is not None
    assert quantifier_instance.quantify is not None
    assert quantifier_instance.calib_weight is not None
    assert quantifier_instance.case in Quantifier.CASES


