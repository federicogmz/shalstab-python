"""
Model Context Protocol (MCP) server for SHALSTAB slope stability analysis.

Exposes deterministic slope stability calculations and hydrological steady-state
diagnostic tools compatible with Claude Desktop, Claude Code, and OpenAI Agents SDK.
"""

from __future__ import annotations

import math
import sys
from typing import Any, Dict

# Water unit weight constant matching Analyzer.GRAVITY_WATER
GRAVITY_WATER: float = 9.81  # kN/m³

# FastMCP / MCPServer import handled gracefully across mcp 1.x and 2.x
try:
    from mcp.server.fastmcp import FastMCP
except ImportError:
    try:
        from mcp.server.mcpserver import MCPServer as FastMCP
    except ImportError:
        FastMCP = None


def compute_infinite_slope_stability(
    slope_degrees: float,
    cohesion_kpa: float,
    friction_angle_degrees: float,
    soil_unit_weight_kn_m3: float = 19.0,
    soil_thickness_m: float = 1.5,
    water_table_ratio: float = 0.0,
) -> Dict[str, Any]:
    """
    Compute factor of safety and SHALSTAB stability classification.

    Evaluates Coulomb limit equilibrium combined with Montgomery & Dietrich (1994)
    infinite slope mechanics across dry, saturated, and steady-state seepage.

    Note on stability codes:
    - Code 1 (Unconditionally Stable): tan(theta) < stable threshold. The slope
      cannot fail even under complete saturation (m=1.0).
    - Code 2 (Unconditionally Unstable): tan(theta) >= dry strength threshold.
      The slope is unstable even in dry conditions (m=0.0).
    - Code 3 (Unstable): conditionally unstable under current water table ratio m.
    - Code 4 (Stable): conditionally stable under current water table ratio m.
    For unconditional categories (codes 1 and 2), the classification is an intrinsic
    property of terrain slope and soil strength, independent of water table height.

    Parameters
    ----------
    slope_degrees : float
        Topographic slope angle theta in degrees (0 < slope_degrees < 90).
    cohesion_kpa : float
        Effective soil cohesion c' in kPa (cohesion_kpa >= 0).
    friction_angle_degrees : float
        Internal friction angle phi' in degrees (0 < friction_angle_degrees < 90).
    soil_unit_weight_kn_m3 : float, default 19.0
        Total soil unit weight gamma in kN/m³ (must be > 9.81 kN/m³).
    soil_thickness_m : float, default 1.5
        Vertical soil mantle thickness z in meters (soil_thickness_m > 0).
    water_table_ratio : float, default 0.0
        Water table height ratio m = h/z (0.0 = dry, 1.0 = fully saturated).

    Returns
    -------
    Dict[str, Any]
        Dictionary containing factor_of_safety, stability_category, stability_code,
        unconditional thresholds, physical diagnosis, and input parameters.
    """
    if not (0.0 < slope_degrees < 90.0):
        raise ValueError(f"slope_degrees must be in (0, 90), got {slope_degrees}")
    if cohesion_kpa < 0.0:
        raise ValueError(f"cohesion_kpa must be >= 0, got {cohesion_kpa}")
    if not (0.0 < friction_angle_degrees < 90.0):
        raise ValueError(
            f"friction_angle_degrees must be in (0, 90), got {friction_angle_degrees}"
        )
    if soil_unit_weight_kn_m3 <= GRAVITY_WATER:
        raise ValueError(
            f"soil_unit_weight_kn_m3 must be > water unit weight ({GRAVITY_WATER}), "
            f"got {soil_unit_weight_kn_m3}"
        )
    if soil_thickness_m <= 0.0:
        raise ValueError(f"soil_thickness_m must be > 0, got {soil_thickness_m}")
    if not (0.0 <= water_table_ratio <= 1.0):
        raise ValueError(
            f"water_table_ratio must be in [0, 1], got {water_table_ratio}"
        )

    theta = math.radians(slope_degrees)
    phi = math.radians(friction_angle_degrees)
    tan_slope = math.tan(theta)
    tan_friction = math.tan(phi)
    cos_slope = math.cos(theta)
    sin_slope = math.sin(theta)
    cos_slope_sq = cos_slope**2

    # Cohesion resistance component
    cohesion_term = cohesion_kpa / (
        soil_unit_weight_kn_m3 * soil_thickness_m * cos_slope_sq
    )

    # Unconditional stability thresholds matching Analyzer
    unstable_threshold = tan_friction + cohesion_term
    stable_threshold = (
        (1.0 - (GRAVITY_WATER / soil_unit_weight_kn_m3)) * tan_friction
    ) + cohesion_term

    # Factor of Safety calculation (Skempton-Delory / Coulomb equation)
    shear_strength = cohesion_kpa + (
        (soil_unit_weight_kn_m3 - water_table_ratio * GRAVITY_WATER)
        * soil_thickness_m
        * cos_slope_sq
        * tan_friction
    )
    shear_stress = soil_unit_weight_kn_m3 * soil_thickness_m * sin_slope * cos_slope
    fs = shear_strength / shear_stress if shear_stress > 0.0 else float("inf")

    # SHALSTAB Classification
    is_unconditionally_unstable = tan_slope >= unstable_threshold
    is_unconditionally_stable = tan_slope < stable_threshold

    if is_unconditionally_unstable:
        category = "Unconditionally Unstable"
        code = 2
        dry_deg = math.degrees(math.atan(unstable_threshold))
        diagnosis = (
            f"Slope angle ({slope_degrees:.1f}°) exceeds dry failure threshold "
            f"({dry_deg:.1f}°). Unstable even in dry conditions (FS = {fs:.2f})."
        )
    elif is_unconditionally_stable:
        category = "Unconditionally Stable"
        code = 1
        sat_deg = math.degrees(math.atan(stable_threshold))
        diagnosis = (
            f"Slope angle ({slope_degrees:.1f}°) is gentler than saturated threshold "
            f"({sat_deg:.1f}°). Stable even under full saturation (FS = {fs:.2f})."
        )
    else:
        if fs < 1.0:
            category = "Unstable"
            code = 3
            diagnosis = (
                f"Conditionally unstable at water table ratio "
                f"m={water_table_ratio:.2f} (FS = {fs:.2f} < 1.0)."
            )
        else:
            category = "Stable"
            code = 4
            diagnosis = (
                f"Conditionally stable at water table ratio "
                f"m={water_table_ratio:.2f} (FS = {fs:.2f} >= 1.0)."
            )

    return {
        "factor_of_safety": round(fs, 4),
        "stability_category": category,
        "stability_code": code,
        "unconditional_stable_threshold_tan": round(stable_threshold, 4),
        "unconditional_unstable_threshold_tan": round(unstable_threshold, 4),
        "is_unconditionally_stable": is_unconditionally_stable,
        "is_unconditionally_unstable": is_unconditionally_unstable,
        "physical_diagnosis": diagnosis,
        "parameters": {
            "slope_degrees": slope_degrees,
            "cohesion_kpa": cohesion_kpa,
            "friction_angle_degrees": friction_angle_degrees,
            "soil_unit_weight_kn_m3": soil_unit_weight_kn_m3,
            "soil_thickness_m": soil_thickness_m,
            "water_table_ratio": water_table_ratio,
        },
    }


def diagnose_topographic_wetness(
    specific_catchment_area_m: float,
    slope_degrees: float,
    transmissivity_m2_day: float,
    rainfall_mm_day: float,
    cohesion_kpa: float = 0.0,
    friction_angle_degrees: float = 33.0,
    soil_unit_weight_kn_m3: float = 19.0,
    soil_thickness_m: float = 1.5,
) -> Dict[str, Any]:
    """
    Diagnose topographic wetness index and critical rainfall failure threshold.

    Computes steady-state hydrologic saturation w = (q * a) / (T * sin(theta))
    and critical log(q/T) ratio per Montgomery & Dietrich (1994).

    Parameters
    ----------
    specific_catchment_area_m : float
        Specific contributing drainage area a = A / b in meters, where A is upslope
        contributing area (m²) and b is contour width (m).
    slope_degrees : float
        Ground surface slope angle theta in degrees (0 < slope_degrees < 90).
    transmissivity_m2_day : float
        Soil profile transmissivity T = K * z * cos(theta) in m²/day.
    rainfall_mm_day : float
        Steady-state rainfall infiltration rate in mm/day. Converted internally
        to recharge q in m/day (q = rainfall_mm_day / 1000.0).
    cohesion_kpa : float, default 0.0
        Effective soil cohesion c' in kPa (cohesion_kpa >= 0.0).
    friction_angle_degrees : float, default 33.0
        Internal friction angle phi' in degrees (0 < friction_angle_degrees < 90).
    soil_unit_weight_kn_m3 : float, default 19.0
        Total soil unit weight gamma in kN/m³ (must be > 9.81 kN/m³).
    soil_thickness_m : float, default 1.5
        Vertical soil depth z in meters (soil_thickness_m > 0).

    Returns
    -------
    Dict[str, Any]
        Hydrologic diagnosis including wetness ratio w, critical recharge ratios
        ('log10_q_over_t_critical' per Montgomery & Dietrich 1994,
        and 'ln_q_over_t_critical' matching Analyzer.calculate_log_qt), failure prediction.
    """
    if specific_catchment_area_m <= 0.0:
        raise ValueError(
            f"specific_catchment_area_m must be > 0, got {specific_catchment_area_m}"
        )
    if not (0.0 < slope_degrees < 90.0):
        raise ValueError(f"slope_degrees must be in (0, 90), got {slope_degrees}")
    if transmissivity_m2_day <= 0.0:
        raise ValueError(
            f"transmissivity_m2_day must be > 0, got {transmissivity_m2_day}"
        )
    if rainfall_mm_day < 0.0:
        raise ValueError(f"rainfall_mm_day must be >= 0, got {rainfall_mm_day}")
    if cohesion_kpa < 0.0:
        raise ValueError(f"cohesion_kpa must be >= 0, got {cohesion_kpa}")
    if not (0.0 < friction_angle_degrees < 90.0):
        raise ValueError(
            f"friction_angle_degrees must be in (0, 90), got {friction_angle_degrees}"
        )
    if soil_unit_weight_kn_m3 <= GRAVITY_WATER:
        raise ValueError(
            f"soil_unit_weight_kn_m3 must be > {GRAVITY_WATER}, "
            f"got {soil_unit_weight_kn_m3}"
        )
    if soil_thickness_m <= 0.0:
        raise ValueError(f"soil_thickness_m must be > 0, got {soil_thickness_m}")

    theta = math.radians(slope_degrees)
    phi = math.radians(friction_angle_degrees)
    sin_slope = math.sin(theta)
    cos_slope = math.cos(theta)
    tan_slope = math.tan(theta)
    tan_friction = math.tan(phi)
    cos_slope_sq = cos_slope**2

    # Recharge q in m/day
    q_m_day = rainfall_mm_day / 1000.0

    # Steady-state wetness ratio w = (q * a) / (T * sin(theta))
    w = (q_m_day * specific_catchment_area_m) / (transmissivity_m2_day * sin_slope)
    w_effective = min(w, 1.0)

    # Check unconditional stability thresholds matching Analyzer
    cohesion_term = cohesion_kpa / (
        soil_unit_weight_kn_m3 * soil_thickness_m * cos_slope_sq
    )
    unstable_threshold = tan_friction + cohesion_term
    stable_threshold = (
        (1.0 - (GRAVITY_WATER / soil_unit_weight_kn_m3)) * tan_friction
    ) + cohesion_term

    is_unconditionally_unstable = tan_slope >= unstable_threshold
    is_unconditionally_stable = tan_slope < stable_threshold

    # Critical recharge threshold (q/T)_crit (Montgomery & Dietrich, 1994)
    term_friction = (soil_unit_weight_kn_m3 / GRAVITY_WATER) * (
        1.0 - (tan_slope / tan_friction)
    )
    term_cohesion = cohesion_kpa / (
        GRAVITY_WATER * soil_thickness_m * cos_slope_sq * tan_friction
    )
    q_over_t_crit = (sin_slope / specific_catchment_area_m) * (
        term_friction + term_cohesion
    )

    if is_unconditionally_unstable:
        log10_q_over_t_crit = None
        ln_q_over_t_crit = None
        failure_predicted = True
        diagnosis = (
            "Slope is unconditionally unstable (tan(theta) exceeds dry strength). "
            "Failure is predicted regardless of hydrological recharge."
        )
    elif is_unconditionally_stable:
        log10_q_over_t_crit = None
        ln_q_over_t_crit = None
        failure_predicted = False
        diagnosis = (
            "Slope is unconditionally stable (tan(theta) is below threshold). "
            "No rainfall amount can trigger shallow planar failure."
        )
    else:
        log10_q_over_t_crit = round(math.log10(q_over_t_crit), 4)
        ln_q_over_t_crit = round(math.log(q_over_t_crit), 4)
        q_over_t_current = q_m_day / transmissivity_m2_day
        failure_predicted = q_over_t_current >= q_over_t_crit
        if failure_predicted:
            diagnosis = (
                f"Current recharge ratio q/T ({q_over_t_current:.6f}) exceeds "
                f"critical threshold ({q_over_t_crit:.6f}). Failure predicted."
            )
        else:
            diagnosis = (
                f"Current recharge ratio q/T ({q_over_t_current:.6f}) is below "
                f"critical threshold ({q_over_t_crit:.6f}). Stable under rainfall."
            )

    if w >= 1.0:
        regime = "Overland flow / Full saturation (w >= 1.0)"
    elif w > 0.5:
        regime = "High pore pressure build-up (0.5 < w < 1.0)"
    else:
        regime = "Subsurface unsaturated regime (w <= 0.5)"

    return {
        "wetness_ratio_w": round(w, 4),
        "effective_saturation_ratio": round(w_effective, 4),
        "saturation_percentage": round(w_effective * 100.0, 2),
        "critical_recharge_ratio_q_over_t": (
            round(q_over_t_crit, 6)
            if (not is_unconditionally_unstable and not is_unconditionally_stable)
            else None
        ),
        "log10_q_over_t_critical": log10_q_over_t_crit,
        "ln_q_over_t_critical": ln_q_over_t_crit,
        "hydrologic_regime": regime,
        "is_saturated": w >= 1.0,
        "is_unconditionally_stable": is_unconditionally_stable,
        "is_unconditionally_unstable": is_unconditionally_unstable,
        "failure_predicted": failure_predicted,
        "physical_diagnosis": diagnosis,
        "parameters": {
            "specific_catchment_area_m": specific_catchment_area_m,
            "slope_degrees": slope_degrees,
            "transmissivity_m2_day": transmissivity_m2_day,
            "rainfall_mm_day": rainfall_mm_day,
            "cohesion_kpa": cohesion_kpa,
            "friction_angle_degrees": friction_angle_degrees,
            "soil_unit_weight_kn_m3": soil_unit_weight_kn_m3,
            "soil_thickness_m": soil_thickness_m,
        },
    }


def create_mcp_server() -> Any:
    """Create and configure the FastMCP server instance."""
    if FastMCP is None:
        raise ImportError(
            "The 'mcp' package is required to run the MCP server. "
            "Install it with: pip install 'shalstab[mcp]'"
        )

    server = FastMCP(
        "shalstab",
        instructions=(
            "Deterministic infinite slope stability and hydrological analysis "
            "based on SHALSTAB (Montgomery & Dietrich, 1994). Evaluates "
            "Factor of Safety and critical rainfall thresholds. "
            "Requires Python >= 3.10."
        ),
    )

    server.tool()(compute_infinite_slope_stability)
    server.tool()(diagnose_topographic_wetness)

    return server


def main() -> None:
    """Entry point for shalstab-mcp CLI command."""
    if FastMCP is None:
        sys.stderr.write(
            "Error: Model Context Protocol SDK ('mcp') is not installed "
            "or current Python version is < 3.10.\n"
            "The SHALSTAB MCP server requires Python >= 3.10. Install with:\n"
            "    pip install 'shalstab[mcp]'\n"
        )
        sys.exit(1)

    server = create_mcp_server()
    server.run()


if __name__ == "__main__":
    main()
