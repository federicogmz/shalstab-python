"""Tests for SHALSTAB Model Context Protocol (MCP) server analytical tools."""

import math

import pytest

from shalstab.mcp_server import (
    compute_infinite_slope_stability,
    diagnose_topographic_wetness,
)


def test_dry_slope_infinite_factor_of_safety():
    """Verify dry slope factor of safety matches tan(phi)/tan(theta)."""
    slope = 20.0
    phi = 35.0
    result = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=0.0,
        friction_angle_degrees=phi,
        water_table_ratio=0.0,
    )
    expected_fs = math.tan(math.radians(phi)) / math.tan(math.radians(slope))
    assert math.isclose(result["factor_of_safety"], expected_fs, rel_tol=1e-3)
    assert result["is_unconditionally_stable"] is False
    assert result["is_unconditionally_unstable"] is False
    assert result["stability_category"] == "Stable"
    assert result["stability_code"] == 4


def test_fully_saturated_factor_of_safety():
    """Verify saturated slope factor of safety matches analytical formula."""
    slope = 15.0
    phi = 32.0
    gamma = 20.0
    gamma_w = 9.81
    result = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=0.0,
        friction_angle_degrees=phi,
        soil_unit_weight_kn_m3=gamma,
        water_table_ratio=1.0,
    )
    expected_fs = (1.0 - (gamma_w / gamma)) * (
        math.tan(math.radians(phi)) / math.tan(math.radians(slope))
    )
    assert math.isclose(result["factor_of_safety"], expected_fs, rel_tol=1e-3)


def test_partially_saturated_intermediate_factor_of_safety():
    """Verify FS decreases monotonically and matches analytical equation at m=0.5."""
    slope = 18.0
    phi = 30.0
    gamma = 19.0
    cohesion = 2.0
    thickness = 1.5
    gamma_w = 9.81

    res_dry = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=cohesion,
        friction_angle_degrees=phi,
        soil_unit_weight_kn_m3=gamma,
        soil_thickness_m=thickness,
        water_table_ratio=0.0,
    )
    res_mid = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=cohesion,
        friction_angle_degrees=phi,
        soil_unit_weight_kn_m3=gamma,
        soil_thickness_m=thickness,
        water_table_ratio=0.5,
    )
    res_sat = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=cohesion,
        friction_angle_degrees=phi,
        soil_unit_weight_kn_m3=gamma,
        soil_thickness_m=thickness,
        water_table_ratio=1.0,
    )

    assert (
        res_dry["factor_of_safety"]
        > res_mid["factor_of_safety"]
        > res_sat["factor_of_safety"]
    )

    # Exact analytical verification for m=0.5
    theta = math.radians(slope)
    phi_rad = math.radians(phi)
    cos_sq = math.cos(theta) ** 2
    shear_strength = cohesion + (
        (gamma - 0.5 * gamma_w) * thickness * cos_sq * math.tan(phi_rad)
    )
    shear_stress = gamma * thickness * math.sin(theta) * math.cos(theta)
    expected_mid_fs = shear_strength / shear_stress
    assert math.isclose(res_mid["factor_of_safety"], expected_mid_fs, rel_tol=1e-3)


def test_unconditionally_unstable_slope():
    """Verify slope steeper than friction angle with zero cohesion is unstable."""
    result = compute_infinite_slope_stability(
        slope_degrees=40.0,
        cohesion_kpa=0.0,
        friction_angle_degrees=30.0,
    )
    assert result["is_unconditionally_unstable"] is True
    assert result["stability_category"] == "Unconditionally Unstable"
    assert result["stability_code"] == 2


def test_unconditionally_stable_slope():
    """Verify gentle slope with high cohesion is stable even when saturated."""
    result = compute_infinite_slope_stability(
        slope_degrees=5.0,
        cohesion_kpa=25.0,
        friction_angle_degrees=35.0,
        water_table_ratio=1.0,
    )
    assert result["is_unconditionally_stable"] is True
    assert result["stability_category"] == "Unconditionally Stable"
    assert result["stability_code"] == 1


def test_invalid_parameters_raise_value_error():
    """Verify physical bounds validation and error messages."""
    with pytest.raises(ValueError, match="slope_degrees"):
        compute_infinite_slope_stability(
            slope_degrees=0.0, cohesion_kpa=0, friction_angle_degrees=30
        )
    with pytest.raises(ValueError, match="cohesion_kpa"):
        compute_infinite_slope_stability(
            slope_degrees=30.0, cohesion_kpa=-1, friction_angle_degrees=30
        )
    with pytest.raises(ValueError, match="friction_angle_degrees"):
        invalid_phi = 95.0
        compute_infinite_slope_stability(
            slope_degrees=30.0, cohesion_kpa=0, friction_angle_degrees=invalid_phi
        )
    with pytest.raises(ValueError, match="soil_unit_weight_kn_m3"):
        compute_infinite_slope_stability(
            slope_degrees=30.0,
            cohesion_kpa=0,
            friction_angle_degrees=30,
            soil_unit_weight_kn_m3=8.0,
        )
    with pytest.raises(ValueError, match="water_table_ratio"):
        compute_infinite_slope_stability(
            slope_degrees=30.0,
            cohesion_kpa=0,
            friction_angle_degrees=30,
            water_table_ratio=1.5,
        )


def test_topographic_wetness_diagnosis():
    """Verify steady-state hydrological diagnosis and critical threshold."""
    result = diagnose_topographic_wetness(
        specific_catchment_area_m=25.0,
        slope_degrees=25.0,
        transmissivity_m2_day=10.0,
        rainfall_mm_day=50.0,
        cohesion_kpa=2.0,
        friction_angle_degrees=32.0,
    )
    assert result["wetness_ratio_w"] > 0
    assert 0 <= result["effective_saturation_ratio"] <= 1.0
    assert result["log10_q_over_t_critical"] is not None
    assert result["ln_q_over_t_critical"] is not None
    # Rigorous mathematical cross-validation between base-10 and natural log
    assert math.isclose(
        result["ln_q_over_t_critical"],
        result["log10_q_over_t_critical"] * math.log(10),
        rel_tol=1e-3,
    )
    assert "cohesion_kpa" in result["parameters"]
    assert "failure_predicted" in result


def test_conditional_slope_rain_induced_failure():
    """Verify conditional slope predicts failure when rainfall exceeds threshold."""
    result = diagnose_topographic_wetness(
        specific_catchment_area_m=50.0,
        slope_degrees=25.0,
        transmissivity_m2_day=5.0,
        rainfall_mm_day=150.0,
        cohesion_kpa=1.0,
        friction_angle_degrees=32.0,
    )
    assert result["failure_predicted"] is True
    assert result["log10_q_over_t_critical"] is not None
    assert result["ln_q_over_t_critical"] is not None


def test_topographic_wetness_unconditional_instability():
    """Verify diagnostic when steep slope fails unconditionally without rain."""
    result = diagnose_topographic_wetness(
        specific_catchment_area_m=10.0,
        slope_degrees=45.0,
        transmissivity_m2_day=10.0,
        rainfall_mm_day=0.0,
        cohesion_kpa=0.0,
        friction_angle_degrees=30.0,
    )
    assert result["is_unconditionally_unstable"] is True
    assert result["critical_recharge_ratio_q_over_t"] is None
    assert result["log10_q_over_t_critical"] is None
    assert result["ln_q_over_t_critical"] is None
    assert result["failure_predicted"] is True


def test_topographic_wetness_unconditional_stability():
    """Verify gentle slope cannot fail regardless of extreme rainfall rate."""
    result = diagnose_topographic_wetness(
        specific_catchment_area_m=1000.0,
        slope_degrees=5.0,
        transmissivity_m2_day=1.0,
        rainfall_mm_day=200.0,
        cohesion_kpa=0.0,
        friction_angle_degrees=33.0,
    )
    assert result["is_unconditionally_stable"] is True
    assert result["failure_predicted"] is False
    assert result["log10_q_over_t_critical"] is None
    assert result["ln_q_over_t_critical"] is None
    assert "unconditionally stable" in result["physical_diagnosis"].lower()


def test_topographic_wetness_validation_errors():
    """Verify input validation on hydrological parameters."""
    with pytest.raises(ValueError, match="specific_catchment_area_m"):
        diagnose_topographic_wetness(
            specific_catchment_area_m=-5.0,
            slope_degrees=20.0,
            transmissivity_m2_day=10.0,
            rainfall_mm_day=10.0,
        )
    with pytest.raises(ValueError, match="transmissivity_m2_day"):
        diagnose_topographic_wetness(
            specific_catchment_area_m=10.0,
            slope_degrees=20.0,
            transmissivity_m2_day=-1.0,
            rainfall_mm_day=10.0,
        )


def test_mcp_server_factory():
    """Verify FastMCP server instance creates tools properly if mcp is installed."""
    pytest.importorskip("mcp")
    import asyncio

    from shalstab.mcp_server import create_mcp_server

    server = create_mcp_server()
    assert server.name == "shalstab"
    tools = asyncio.run(server.list_tools())
    tool_names = [t.name for t in tools]
    assert "compute_infinite_slope_stability" in tool_names
    assert "diagnose_topographic_wetness" in tool_names


def test_unconditional_thresholds_match_analytical_formulation():
    """Verify unconditional thresholds match Montgomery & Dietrich (1994) equations."""
    slope = 28.0
    c = 3.5
    phi = 34.0
    gamma = 19.5
    z = 1.2
    gw = 9.81

    res = compute_infinite_slope_stability(
        slope_degrees=slope,
        cohesion_kpa=c,
        friction_angle_degrees=phi,
        soil_unit_weight_kn_m3=gamma,
        soil_thickness_m=z,
    )

    theta = math.radians(slope)
    phi_rad = math.radians(phi)
    cos_sq = math.cos(theta) ** 2
    c_term = c / (gamma * z * cos_sq)
    expected_unstable = math.tan(phi_rad) + c_term
    expected_stable = ((1.0 - (gw / gamma)) * math.tan(phi_rad)) + c_term

    assert math.isclose(
        res["unconditional_unstable_threshold_tan"], expected_unstable, rel_tol=1e-4
    )
    assert math.isclose(
        res["unconditional_stable_threshold_tan"], expected_stable, rel_tol=1e-4
    )
