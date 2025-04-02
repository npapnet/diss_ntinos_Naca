import pytest

import algorithmos_DTU_10_sections as d10

@pytest.fixture
def bl_cl():
    blade_geom_DTU = "blade_geom_DTU.json"
    hansen_DTU = d10.Hansen_Algorithm_for_DTU_geometry(
        blade_geom_DTU=blade_geom_DTU,
        B=3,
        air_density=1.225
    )
    return hansen_DTU

def test_object_creation(bl_cl):
    assert bl_cl.B == 3
    assert bl_cl.R  == pytest.approx(89.17, abs=0.010)

def test_calculation_of_flow_angle(bl_cl):
    v0=10
    omega_rps=0.5
    actual = bl_cl.calculation_of_flow_angle_rad(a=0, a_p=0, r=2.8, v0=v0, w_rps=omega_rps)
    expected = 1.431 
    assert actual == pytest.approx(expected, rel=1e-2)
    
def test_calculation_of_local_angle_of_attack(bl_cl):
    pitch_angle_deg = 14.5
    twist_deg = 0 
    actual = bl_cl.calculation_of_local_angle_of_attack_rad(flow_angle_rad=1.431, pitch_angle_deg=pitch_angle_deg, twist_deg=twist_deg)
    expected = 1.178 
    assert actual == pytest.approx(expected, rel=1e-2)
    
def test_calculation_of_Cl_and_Cd(bl_cl):
    angle_of_attack_deg = 67.53
    tc_ratio = 100
    actual = bl_cl.calcualtion_of_Cl_and_Cd(angle_of_attack_deg=angle_of_attack_deg, tc_ratio=tc_ratio)
    expected_Cl, expected_Cd = 0, 0.6
    assert actual == pytest.approx(expected_Cd, expected_Cd, rel=1e-2)
    
def test_calcualtion_of_Cn_and_Ct(bl_cl):
    Cl = 0
    Cd = 0.6
    flow_angle_rad = 1.431
    actual = bl_cl.calculation_of_Cn_and_Ct(Cl=Cl, Cd=Cd, flow_angle_rad=flow_angle_rad)
    expected_Cn, expected_Ct = 0.594, -0.08
    assert actual == pytest.approx(expected_Cn, expected_Ct, rel=1e-2)
    