import pytest

import algorithmos_DTU_10_sections as d10

@pytest.fixture
def bl_cl():
    blade_geom_file_2 = "blade_geom_DTU.json"
    hansen_DTU = d10.Hansen_Algorithm_for_DTU_geometry(
        blade_geom_file_2=blade_geom_file_2,
        B=3,
        air_density=1.225
    )
    return hansen_DTU

def test_object_creation(bl_cl):
    assert bl_cl.B == 3
    assert bl_cl.R  == pytest.approx(89.17, abs=0.010)
    

def test_calculate_flow_angle(bl_cl):
    v0=10
    omega_rps=0.5
    actual = bl_cl.calculation_of_flow_angle_rad(a=0,ap=0, r=2.8,v0=v0, w_rps=omega_rps )
    expected = 0.12 
    assert actual == pytest.approx(expected, rel=1e-2)
    
    
    