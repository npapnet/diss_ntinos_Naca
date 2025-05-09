import pytest

import _algorithmos_DTU as d10

@pytest.fixture
def bl_cl():
    blade_geom_DTU = "blade_geom_DTU.json"
    hansen_DTU = d10.Hansen_Algorithm(
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
    flow_angle_rad = 1.431
    pitch_angle_deg = 14.5
    twist_deg = 0 
    actual = bl_cl.calculation_of_local_angle_of_attack_rad(flow_angle_rad=flow_angle_rad, pitch_angle_deg=pitch_angle_deg, twist_deg=twist_deg)
    expected = 1.178 
    assert actual == pytest.approx(expected, rel=1e-2)
    
def test_calculation_of_Cl_and_Cd(bl_cl):
    angle_of_attack_deg = 67.53
    tc_ratio = 100
    Cl_actual, Cd_actual = bl_cl.calculation_of_Cl_and_Cd(angle_of_attack_deg=angle_of_attack_deg, tc_ratio=tc_ratio)
    expected_Cl = 0
    expected_Cd = 0.6
    assert Cl_actual == pytest.approx(expected_Cl, rel=1e-2)
    assert Cd_actual == pytest.approx(expected_Cd, rel=1e-2)
    
def test_calcualtion_of_Cn_and_Ct(bl_cl):
    Cl = 0
    Cd = 0.6
    flow_angle_rad = 1.431
    Cn_actual, Ct_actual = bl_cl.calculation_of_Cn_and_Ct(Cl, Cd, flow_angle_rad)
    expected_Cn = 0.594
    expected_Ct = -0.08
    assert Cn_actual == pytest.approx(expected_Cn, rel=1e-1)
    assert Ct_actual == pytest.approx(expected_Ct, rel=1e-1)
    
def test_calculation_of_updated_induction_factors(bl_cl):
    Cn = 0.594
    Ct = -0.08
    r = 2.8
    chord = 5.38
    flow_angle_rad = 1.431
    
    a_new, a_p_new = bl_cl.calculation_of_updated_induction_factors(Cn=Cn, Ct=Ct, r=r, chord=chord, flow_angle_rad=flow_angle_rad)
    expected_a_new = 0.122
    expected_a_p_new = -0.122
    assert a_new == pytest.approx(expected_a_new, rel=1e-1)
    assert a_p_new == pytest.approx(expected_a_p_new, rel=1e-1)
    
def test_calculation_of_local_loads(bl_cl):
    r = 2.8
    a, a_p = 0, 0
    v0 = 10.0
    w_rps = 0.5
    chord = 5.38
    flow_angle_rad = 1.431
    Cl, Cd = 0, 0.6
    
    L_actual, D_actual, pn_actual, pt_actual = bl_cl.calculation_of_local_loads(
        r=r, a=a, a_p=a_p, v0=v0, w_rps=w_rps, chord=chord, flow_angle_rad=flow_angle_rad, Cl=Cl, Cd=Cd)
    
    expected_L = 0
    expected_D = 201.59
    expected_pn = 199.64
    expected_pt = -27.95
    
    assert L_actual == pytest.approx(expected_L, rel=1e-2)
    assert D_actual == pytest.approx(expected_D, rel=1e-2)
    assert pn_actual == pytest.approx(expected_pn, rel=1e-2)
    assert pt_actual == pytest.approx(expected_pt, rel=1e-2)

    
    