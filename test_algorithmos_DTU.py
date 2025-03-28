import pytest

import algorithmos_DTU_10_sections as d10


def test_one():
    blade_geom_file_2 = "blade_geom_file_2.json"
    hansen_DTU = d10.Hansen_Algorithm_for_DTU_geometry(
        blade_geom_file_2=blade_geom_file_2,
        B=3,
        air_density=1.225
    )
    