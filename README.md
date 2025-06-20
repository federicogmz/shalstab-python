# SHALSTAB - Shallow Landsliding STABility Model

A comprehensive Python package for slope stability analysis using the SHALSTAB (Shallow Landsliding STABility) model. This implementation provides tools for analyzing slope stability using physically-based models that consider topographic, hydrologic, and geotechnical factors.

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Scientific Background](#scientific-background)
- [Mathematical Framework](#mathematical-framework)
- [Model Components](#model-components)
- [Usage Examples](#usage-examples)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

## Overview

The SHALSTAB model evaluates infinite slope stability by combining:

- **Topographic analysis** (slope, flow accumulation)
- **Hydrologic modeling** (steady-state water flow)
- **Geotechnical parameters** (cohesion, friction angle, unit weight, permeability)

### Key Features

- **Complete SHALSTAB implementation** with all core analyses
- **Raster-based processing** using xarray and rasterio
- **Professional visualizations** with matplotlib
- **Comprehensive error handling** and validation
- **Type hints** and detailed documentation
- **Export capabilities** to GeoTIFF and other formats

## Installation

### Using pip (when available)

```bash
pip install shalstab
```

### From source

```bash
git clone https://github.com/federicogmz/shalstab.git
cd shalstab
pip install -e .
```

### Dependencies

```bash
pip install -r requirements.txt
```

## Training data

The easiest way to test if SHALSTAB is installed correctly or to get familiar with its implementation is to use the provided training data:

```python
import shalstab

# Use training data
analyzer = shalstab.Analyzer(
    dem_path=shalstab.training_dem,
    geo=shalstab.training_geology,
    geo_columns=["Cohesion", "Phi", "Gamma_kN_m", "Ks_m_s"]
)

# Run analysis
critical_rain = analyzer.calculate_critical_rainfall()
```

## Scientific Background

### SHALSTAB Model Overview

SHALSTAB (Shallow Landsliding STABility) is a physically-based model developed by Montgomery and Dietrich (1994) for analyzing shallow landslide susceptibility. The model combines:

1. **Infinite slope stability analysis** - Classical geotechnical stability assessment
2. **Steady-state hydrologic model** - Water flow and saturation analysis
3. **Topographic controls** - Digital terrain analysis for slope and drainage

The model is particularly suited for analyzing rainfall-triggered shallow landslides in natural slopes where the failure surface is parallel to the ground surface.

### Model Assumptions

The SHALSTAB model is based on several key assumptions:

1. **Infinite slope geometry**: The slope extends infinitely in all directions, allowing for 1D analysis
2. **Planar failure surface**: Failure occurs along a plane parallel to the ground surface
3. **Steady-state hydrology**: Groundwater flow reaches equilibrium during rainfall events
4. **Saturated conditions**: The model considers fully saturated soil conditions at failure
5. **Uniform soil properties**: Geotechnical parameters are assumed constant within each geological unit
6. **Coulomb failure criterion**: Soil strength follows the Mohr-Coulomb failure envelope

### Applications and Limitations

**Suitable for:**

- Regional-scale landslide susceptibility mapping
- Rainfall threshold analysis
- Comparative stability assessment
- Risk prioritization and hazard zoning

**Limitations:**

- Not suitable for deep-seated landslides
- Cannot predict exact timing of failure
- Requires high-quality input data
- May not capture complex hydrogeological conditions

## Mathematical Framework

### Core SHALSTAB Equation

The fundamental SHALSTAB stability criterion compares the ratio of upslope contributing area to contour length (destabilizing force) with a stabilizing term that includes soil properties and rainfall:

```
a/b ≤ (T/q) × [(c/(γw × z × cos²θ × tanφ)) + (γ/γw) × (1 - tanθ/tanφ)]
```

Where:

- `a` = upslope contributing area (m²)
- `b` = contour length (cell width, m)
- `T` = transmissivity (k × z × cosθ, m²/day)
- `q` = recharge rate (rainfall intensity, mm/day converted to m/day)
- `c` = soil cohesion (kPa)
- `γw` = unit weight of water (9.81 kN/m³)
- `z` = soil thickness (m)
- `θ` = slope angle (radians)
- `φ` = internal friction angle (radians)
- `γ` = soil unit weight (kN/m³)
- `k` = saturated hydraulic conductivity (m/s)

### Critical Rainfall Calculation

The critical rainfall represents the minimum rainfall intensity required to trigger slope failure:

```
qcrit = (T/a) × [(c/(γw × z × cos²θ × tanφ)) + (γ/γw) × (1 - tanθ/tanφ)]
```

This equation is derived by rearranging the stability criterion to solve for the critical recharge rate.

### Stability Classification

The model classifies each cell into one of four stability categories:

1. **Unconditionally Stable (Class 1)**: `tanθ/tanφ < γ/γw` - Stable regardless of rainfall
2. **Unconditionally Unstable (Class 2)**: `tanθ/tanφ > 1` - Unstable regardless of rainfall
3. **Unstable (Class 3)**: Left side > Right side for given rainfall - Unstable at specified rainfall
4. **Stable (Class 4)**: Left side ≤ Right side for given rainfall - Stable at specified rainfall

### Log(q/T) Analysis

The log(q/T) ratio provides a dimensionless measure of stability:

```
log(q/T) = log10(q × cell_size / flow_accumulated) - log10(T)
```

Where:

- Negative values indicate potentially unstable conditions
- Positive values indicate more stable conditions
- Values near zero represent marginal stability

## Model Components

### 1. Soil Thickness Estimation

The package implements the Catani et al. (2010) empirical model for soil thickness estimation:

```
z = zmax × [1 - ((tanθ - tanθmin) / (tanθmax - tanθmin)) × (1 - zmin/zmax)]
```

Where:

- `z` = calculated soil thickness (m)
- `zmax` = maximum soil thickness (typically 5.0 m)
- `zmin` = minimum soil thickness (typically 0.1 m)
- `θ` = local slope angle (radians)
- `θmax`, `θmin` = maximum and minimum slope angles in study area

**Model Assumptions:**

- Linear relationship between slope and soil thickness
- Soil loss increases with slope steepness
- Minimum thickness maintained even on steepest slopes
- Maximum thickness occurs on flattest areas

**Applications:**

- Weight of soil column calculation (driving force)
- Available cohesive strength estimation
- Transmissivity calculations (combined with permeability)

### 2. Topographic Analysis

#### Slope Calculation

Slope angles are calculated using finite difference approximation:

```
slope = arctan(√((dz/dx)² + (dz/dy)²))
```

The calculation process involves:

1. Computing elevation gradients using finite differences
2. Calculating slope magnitude from gradient components
3. Converting from radians to degrees for display

#### Flow Direction and Accumulation

The package uses the D8 flow direction algorithm:

1. **Flow Direction**: Determines steepest descent direction for each cell
2. **Flow Accumulation**: Calculates upslope contributing area using D8 routing
3. **Preprocessing**: Fills pits and resolves flat areas for continuous flow routing

### 3. Hydrologic Modeling

#### Transmissivity Calculation

Soil transmissivity represents the ability to transmit water laterally:

```
T = k × z × cosθ
```

Where:

- `T` = transmissivity (m²/day)
- `k` = saturated hydraulic conductivity (m/s, converted to m/day)
- `z` = soil thickness (m)
- `θ` = slope angle (radians)

#### Steady-State Flow Model

The model assumes steady-state conditions where:

- Recharge rate equals rainfall intensity
- Groundwater flow reaches equilibrium
- Saturated zone thickness remains constant during analysis

### 4. Geotechnical Framework

#### Mohr-Coulomb Failure Criterion

Soil strength follows the Mohr-Coulomb relationship:

```
τ = c + σn × tanφ
```

Where:

- `τ` = shear strength (kPa)
- `c` = cohesion (kPa)
- `σn` = normal stress (kPa)
- `φ` = internal friction angle (degrees)

#### Infinite Slope Analysis

For infinite slope conditions, the factor of safety is:

```
FS = (c + (γ×z×cosθ - γw×zw×cosθ) × tanφ) / (γ×z×sinθ×cosθ)
```

Where:

- `FS` = factor of safety
- `zw` = height of water table above failure plane (m)

### 5. Failure Probability Assessment

The failure probability calculation provides a relative ranking of landslide susceptibility:

**Process:**

1. Calculate log(q/T) using `calculate_log_qt()`
2. Invert values: `-log(q/T)`
3. Normalize to 0-100% range: `P = ((inverted - min) / (max - min)) × 100`

**Interpretation:**

- Areas requiring high rainfall to fail → Low probability
- Areas requiring low rainfall to fail → High probability
- Unconditionally stable/unstable areas → No probability (NaN)

**Important Note:** The probability is relative, not absolute. It represents comparative susceptibility across the study area rather than actual failure rates.

## Usage Examples

### Basic Topographic Analysis

```python
import shalstab

# Example 1: Using file paths (recommended)
analyzer = shalstab.Analyzer(
    dem_path="elevation.tif",
    geo="geology.geojson",
    geo_columns=["cohesion", "friction", "gamma", "permeability"]
)

# Example 2: Using GeoDataFrame
import geopandas as gpd
geology = gpd.read_file("geology.shp")
analyzer = shalstab.Analyzer("elevation.tif", geology, geo_columns=["cohesion", "friction", "gamma", "permeability"])
```

```python
# Initialize full analyzer with file paths
analyzer = shalstab.Analyzer(
    dem_path="high_res_dem.tif",
    geo="geology_with_properties.geojson",
    geo_columns=[
        "cohesion_kpa",      # Soil cohesion (kPa)
        "friction_deg",      # Internal friction angle (degrees)
        "gamma_knm3",        # Soil unit weight (kN/m³)
        "k_ms"               # Hydraulic conductivity (m/s)
    ],
    figsize=(15, 10)
)

# Calculate critical rainfall with visualization
critical_rainfall = analyzer.calculate_critical_rainfall(show_plot=True)
print(f"Critical rainfall range: {critical_rainfall.min():.1f} - {critical_rainfall.max():.1f} mm/day")

# Analyze stability for specific rainfall events
rainfall_events = [10, 25, 50, 100]  # mm/day

for rainfall in rainfall_events:
    stability, fig = analyzer.calculate_stability(rainfall_mm_day=rainfall)

    # Calculate unstable area
    unstable_cells = (stability == 3).sum()
    cell_area = abs(analyzer.dem.rio.resolution()[0] * analyzer.dem.rio.resolution()[1])
    unstable_area_km2 = unstable_cells * cell_area / 1e6

    print(f"Rainfall {rainfall} mm/day: {unstable_area_km2:.2f} km² unstable")

    # Save results
    analyzer.export_raster(stability, f"stability_{rainfall}mm.tif")
    fig.savefig(f"stability_plot_{rainfall}mm.png", dpi=300, bbox_inches='tight')

# Calculate relative failure probability
probability = analyzer.calculate_failure_probability()
print(f"High risk areas (>80%): {(probability > 80).sum()} cells")

# Export all results
analyzer.export_raster(critical_rainfall, "critical_rainfall.tif")
analyzer.export_raster(probability, "failure_probability.tif")
analyzer.export_raster(analyzer.soil_thickness, "soil_thickness.tif")
```

### Advanced Analysis Workflow

```python
# Multi-scenario analysis
scenarios = {
    "dry_season": 15,      # mm/day
    "normal_rain": 35,     # mm/day
    "heavy_rain": 75,      # mm/day
    "extreme_event": 150   # mm/day
}

results = {}
for scenario, rainfall in scenarios.items():
    stability, fig = analyzer.calculate_stability(rainfall_mm_day=rainfall)

    # Calculate statistics
    stats = {
        'total_cells': stability.count().item(),
        'stable_cells': (stability == 4).sum().item(),
        'unstable_cells': (stability == 3).sum().item(),
        'uncond_stable': (stability == 1).sum().item(),
        'uncond_unstable': (stability == 2).sum().item()
    }

    # Calculate percentages
    total = stats['total_cells']
    stats['stable_pct'] = stats['stable_cells'] / total * 100
    stats['unstable_pct'] = stats['unstable_cells'] / total * 100

    results[scenario] = stats

    # Export scenario
    analyzer.export_raster(stability, f"stability_{scenario}.tif")
    fig.savefig(f"stability_{scenario}.png", dpi=300, bbox_inches='tight')

# Print summary
print("\nStability Analysis Summary:")
print("-" * 60)
for scenario, stats in results.items():
    print(f"{scenario:15} | Stable: {stats['stable_pct']:5.1f}% | Unstable: {stats['unstable_pct']:5.1f}%")
```

### Data Preprocessing

```python
# Preprocess DEM to fill NoData values
processed_dem = analyzer.preprocess_dem("raw_dem_with_gaps.tif")
print(f"Processed DEM shape: {processed_dem.shape}")

# Resample raster to different resolution
coarse_dem = analyzer.resample_raster("high_res_dem.tif", resolution=30.0)  # 30m resolution
analyzer.export_raster(coarse_dem, "dem_30m.tif")
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.

## License

GNU General Public License v3.0 (GPL-3.0)

See [LICENSE](LICENSE) for full license text.

## References

- Montgomery, D.R., & Dietrich, W.E. (1994). A physically based model for the topographic control on shallow landsliding. Water Resources Research, 30(4), 1153–1171.
- Dietrich, W.E., & Montgomery, D.R. (1998). SHALSTAB: A digital terrain model for mapping shallow landslide potential. NCASI Technical Bulletin 796.
- Catani, F., Segoni, S., & Falorni, G. (2010). An empirical geomorphology-based approach to the spatial prediction of soil thickness at catchment scale. Water Resources Research, 46(5).
