"""
SHALSTAB

This module provides a comprehensive implementation of the SHALSTAB
(Shallow Landsliding STABility) model for slope stability analysis.
"""

import rasterio
import rioxarray
import numpy as np
import xarray as xr
import geopandas as gpd
from pathlib import Path
from scipy import ndimage
from pysheds.grid import Grid
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from geocube.api.core import make_geocube
from typing import Optional, List, Tuple, Union
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap, BoundaryNorm


class Analyzer:
    """
    A comprehensive SHALSTAB slope stability analyzer.

    The SHALSTAB model evaluates infinite slope stability by combining:
    - Topographic analysis (slope, flow accumulation)
    - Hydrologic modeling (steady-state water flow)
    - Geotechnical parameters (cohesion, friction angle, unit weight, permeability)

    Attributes:
        GRAVITY_WATER (float): Water unit weight constant (9.81 kN/m³)
        DEFAULT_FIGSIZE (tuple): Default figure size for plots (12, 8)
        NODATA_VALUE (int): NoData value for raster outputs (-9999)
        STABILITY_CATEGORIES (dict): Mapping of stability codes to descriptive labels
        STABILITY_COLORS (dict): Color scheme for stability visualization
        PROBABILITY_COLORS (list): Color scheme for probability visualization
        CRITICAL_COLORS (list): Color scheme for critical rainfall visualization
    """

    # Class constants
    GRAVITY_WATER = 9.81  # Water unit weight (kN/m³)
    DEFAULT_FIGSIZE = (12, 8)
    NODATA_VALUE = -9999
    STABILITY_CATEGORIES = {
        1: "Unconditionally Stable",
        2: "Unconditionally Unstable",
        3: "Unstable",
        4: "Stable",
    }

    # Color schemes
    STABILITY_COLORS = {1: "darkgreen", 2: "darkred", 3: "gold", 4: "lightgreen"}

    PROBABILITY_COLORS = ["#FFE7F9", "#85E2FF", "#2786EB", "#6A0BA8", "#DE077D"]
    CRITICAL_COLORS = ["#DE077D", "#6A0BA8", "#2786EB", "#85E2FF", "#FFE7F9"]

    def __init__(
        self,
        dem_path: Union[str, Path],
        geo: Union[gpd.GeoDataFrame, str, Path],
        geo_columns: List[str],
        figsize: Tuple[int, int] = DEFAULT_FIGSIZE,
    ):
        """
        Initialize the SHALSTAB analyzer with required datasets and parameters.

        This constructor sets up the analysis environment by loading the Digital Elevation Model,
        processing geological data, and initializing base hydrological and topographic calculations.

        Parameters
        ----------
        dem_path : Union[str, pathlib.Path]
            Path to the Digital Elevation Model raster file. Supported formats include
            GeoTIFF (.tif), ESRI ASCII Grid (.asc), and other GDAL-readable formats.
            The DEM should be in a projected coordinate system with consistent units.

        geo : Union[str, pathlib.Path, geopandas.GeoDataFrame]
            Path to a geological shapefile or GeoJSON, or a GeoDataFrame containing geological units.
            Must include geometry column with polygon features and required attribute columns.

        geo_columns : List[str]
            List of exactly 4 column names in geo GeoDataFrame specifying geotechnical parameters:
            [cohesion_column, friction_angle_column, unit_weight_column, permeability_column]
            - cohesion_column: Soil cohesion values (kPa or kN/m²)
            - friction_angle_column: Internal friction angle (degrees)
            - unit_weight_column: Soil unit weight (kN/m³)
            - permeability_column: Saturated hydraulic conductivity (m/s)

        figsize : Tuple[int, int], default (12, 8)
            Figure size for generated plots as (width, height) in inches.

        Raises
        ------
        ValueError
            If DEM file cannot be loaded or geo_columns has incorrect length.
        FileNotFoundError
            If dem_path does not exist.
        KeyError
            If specified columns in geo_columns are not found in geo GeoDataFrame.

        Attributes Created
        ------------------
        dem : xarray.DataArray
            Loaded Digital Elevation Model with spatial coordinates and CRS information.
        extent : List[float]
            Spatial extent as [min_x, max_x, min_y, max_y] for plotting.
        hillshade : numpy.ndarray
            Calculated hillshade for terrain visualization.
        flow_accumulated : xarray.DataArray
            Flow accumulation grid (m²) calculated using D8 flow direction algorithm.
        slope_rad : xarray.DataArray
            Slope angle in radians calculated from DEM using finite differences.
        cohesion : xarray.DataArray, optional
            Rasterized soil cohesion values (kPa). Only if geo_columns provided.
        friction_rad : xarray.DataArray, optional
            Rasterized friction angle in radians. Only if geo_columns provided.
        unit_weight : xarray.DataArray, optional
            Rasterized soil unit weight (kN/m³). Only if geo_columns provided.
        permeability : xarray.DataArray, optional
            Rasterized hydraulic conductivity (m/s). Only if geo_columns provided.
        _soil_thickness : xarray.DataArray, optional
            Calculated soil thickness using Catani model (m). Only if geo_columns provided.
        unconditional_stability : numpy.ndarray, optional
            Unconditional stability classification. Only if geo_columns provided.
        """
        self.figsize = figsize
        self.dem_path = Path(dem_path)
        # Load geological data: accept file path or GeoDataFrame
        if isinstance(geo, (str, Path)):
            geo_path = Path(geo)
            if not geo_path.exists():
                raise FileNotFoundError(f"Geology file not found: {geo_path}")
            geo_gdf = gpd.read_file(geo_path)
        elif isinstance(geo, gpd.GeoDataFrame):
            geo_gdf = geo
        else:
            raise ValueError(
                "`geo` must be a GeoDataFrame or a path to a shapefile/GeoJSON file"
            )
        self.geo = geo_gdf

        # Enforce required geo_columns argument
        if geo_columns is None:
            raise ValueError("geo_columns is required and cannot be None")

        # Load and validate DEM
        self._load_dem()

        # Calculate extent from geodataframe bounds
        self.extent = self._calculate_extent()

        # Initialize base calculations
        self._initialize_base_calculations()
        # Initialize required geotechnical parameters
        self._initialize_geotechnical_parameters(geo_columns)

    # =====================================================
    # INITIALIZATION AND SETUP METHODS
    # =====================================================

    def _load_dem(self) -> None:
        """Load and validate the Digital Elevation Model."""
        try:
            self.dem = rioxarray.open_rasterio(str(self.dem_path))
        except Exception as e:
            raise ValueError(f"Failed to load DEM from {self.dem_path}: {e}")

    def _calculate_extent(self) -> List[float]:
        """Calculate spatial extent from geodataframe bounds."""
        bounds = self.geo.total_bounds
        return [bounds[0], bounds[2], bounds[1], bounds[3]]

    def _initialize_base_calculations(self) -> None:
        """Initialize base calculations needed for analysis."""
        self.hillshade = self._calculate_hillshade()
        self.flow_accumulated = self._calculate_flow_accumulated()
        self.slope_rad = self._calculate_slope()

    def _initialize_geotechnical_parameters(self, geo_columns: List[str]) -> None:
        """Initialize geotechnical parameters from geodataframe columns."""
        if len(geo_columns) != 4:
            raise ValueError(
                "geo_columns must contain exactly 4 elements: [cohesion, friction, gamma, permeability]"
            )

        cohesion, friction, gamma, permeability = geo_columns
        self.cohesion = self._rasterize_column(cohesion)
        friction_degrees = self._rasterize_column(friction)
        self.friction_rad = np.deg2rad(friction_degrees)
        self.unit_weight = self._rasterize_column(gamma)
        self.permeability = self._rasterize_column(permeability)
        # Calculate soil thickness and stability
        self._soil_thickness = self._calculate_soil_thickness_catani()
        self.unconditional_stability = self._calculate_unconditional_stability()

    # =====================================================
    # DATA PROCESSING AND UTILITY METHODS
    # =====================================================

    def preprocess_dem(self, dem_path: Union[str, Path]) -> xr.DataArray:
        """
        Preprocess Digital Terrain Model by filling NoData values using spatial interpolation.

        This method improves DEM quality by identifying and filling NoData pixels
        (holes) using values from neighboring cells. It applies binary dilation
        to expand valid data regions and uses mean interpolation for gap filling.

        Parameters
        ----------
        dem_path : Union[str, pathlib.Path]
            Path to the Digital Elevation Model raster file to be processed.
            Supported formats include GeoTIFF (.tif), ESRI ASCII Grid (.asc),
            and other GDAL-readable raster formats.

        Returns
        -------
        xarray.DataArray
            Preprocessed DEM with the following properties:
            - NoData values filled using neighboring cell averages
            - Same spatial extent and resolution as input
            - Preserved coordinate reference system
            - Enhanced data continuity for hydrological analysis
            - dtype: float64 (elevation values)
            - Units: Same as input DEM (typically meters)

        Raises
        ------
        FileNotFoundError
            If dem_path does not exist or is not accessible.
        PermissionError
            If the file cannot be overwritten (read-only permissions).
        rasterio.errors.RasterioError
            If the file format is not supported or corrupted.
        """
        print("--- Preprocessing DTM ---")

        dem_path = Path(dem_path)

        with rasterio.open(dem_path) as src:
            data = src.read(1, masked=True)
            nodata_mask = data.mask

            if nodata_mask.any():
                dilated_mask = ndimage.binary_dilation(nodata_mask)
                data_filled = data.copy()
                data_filled[nodata_mask] = np.mean(data[dilated_mask])

                profile = src.profile
                with rasterio.open(dem_path, "w", **profile) as dst:
                    dst.write(data_filled, 1)

        print("--- DTM successfully processed ---")
        return xr.open_dataarray(dem_path)

    def _rasterize_column(self, column_name: str) -> xr.DataArray:
        """
        Convert a vector column to raster format matching the DEM grid.

        This internal method rasterizes attribute data from the geological
        GeoDataFrame onto the same spatial grid as the Digital Elevation Model.
        It uses area-weighted averaging for overlapping polygons and maintains
        spatial consistency across all geotechnical parameter layers.

        Parameters
        ----------
        column_name : str
            Name of the attribute column in self.geo GeoDataFrame to rasterize.
            Must exist in the GeoDataFrame columns.

        Returns
        -------
        xarray.DataArray
            Rasterized attribute data with properties:
            - Same spatial extent and resolution as DEM
            - Same coordinate reference system as DEM
            - dtype: float64 (attribute values)
            - NoData: Areas outside geological polygons filled with NaN
            - Coordinates: x (longitude/easting), y (latitude/northing)

        Raises
        ------
        ValueError
            If column_name does not exist in the GeoDataFrame columns.
        RuntimeError
        If rasterization fails due to memory limitations or incompatible coordinate systems.

        This method is called internally during geotechnical parameter
        initialization and should not typically be called directly by users.
        """
        if column_name not in self.geo.columns:
            raise ValueError(f"Column '{column_name}' not found in geodataframe")

        return make_geocube(
            vector_data=self.geo, measurements=[column_name], like=self.dem, fill=np.nan
        )[column_name]

    def resample_raster(
        self,
        input_path: Union[str, Path],
        resolution: float,
    ) -> xr.DataArray:
        """
        Resample raster to a new resolution.

        Parameters
        ----------
        input_path : Union[str, pathlib.Path]
            Path to input raster file to be resampled.
        resolution : float
            Target resolution for resampling in the same units as the input raster.

        Returns
        -------
        xarray.DataArray
            Resampled raster as DataArray with the specified resolution.
        """
        input_raster = rioxarray.open_rasterio(input_path)

        resampled = input_raster.rio.reproject(
            dst_crs=input_raster.rio.crs, resolution=resolution
        )[0]

        return resampled

    def export_raster(
        self,
        data: xr.DataArray,
        output_path: Union[str, Path],
        driver: str = "GTiff",
        compress: str = "lzw",
        nodata: Optional[float] = None,
        **kwargs,
    ) -> None:
        """
        Export DataArray to raster file using rioxarray.

        Parameters
        ----------
        data : xarray.DataArray
            Input raster data to export with spatial coordinates.
        output_path : Union[str, pathlib.Path]
            Path for the output raster file.
        driver : str, default "GTiff"
            GDAL driver for output format (e.g., "GTiff", "netCDF", "ENVI").
        compress : str, default "lzw"
            Compression method for GeoTIFF files.
        nodata : Optional[float], default None
            NoData value. If None, uses class NODATA_VALUE.
        **kwargs
            Additional keyword arguments passed to rio.to_raster().

        Returns
        -------
        None
            Writes raster file to disk.
        """
        output_path = Path(output_path)

        if nodata is None:
            nodata = self.NODATA_VALUE

        # Ensure data has CRS information
        if data.rio.crs is None and hasattr(self, "dem"):
            data = data.rio.write_crs(self.dem.rio.crs)

        data.rio.to_raster(
            output_path, driver=driver, compress=compress, nodata=nodata, **kwargs
        )

        print(f"Raster exported successfully to {output_path}")

    def _calculate_soil_thickness_catani(
        self, h_min: Union[float, str] = 0.2, h_max: Union[float, str] = 4.5
    ) -> xr.DataArray:
        """
        Calculate soil thickness using Catani et al. (2010) model.

        Args:
            h_min: Minimum soil thickness (constant or column name)
            h_max: Maximum soil thickness (constant or column name)

        Returns:
            Soil thickness as DataArray
        """
        # Handle string inputs (column names)
        if isinstance(h_min, str):
            h_min = self._rasterize_column(h_min)
        if isinstance(h_max, str):
            h_max = self._rasterize_column(h_max)

        # Calculate tangent values
        tan_slope = np.tan(self.slope_rad)
        tan_slope_max = np.tan(np.nanmax(self.slope_rad))
        tan_slope_min = np.tan(np.nanmin(self.slope_rad))

        # Calculate soil thickness
        thickness = h_max * (
            1
            - ((tan_slope - tan_slope_min) / (tan_slope_max - tan_slope_min))
            * (1 - (h_min / h_max))
        )

        thickness = xr.DataArray(
            thickness, coords=[self.dem.coords["y"], self.dem.coords["x"]]
        )

        return thickness.where(~np.isnan(self.dem), np.nan).squeeze()

    def _calculate_slope(self) -> xr.DataArray:
        """
        Calculate slope from DEM.

        Returns:
            Slope in radians as DataArray
        """

        dx, dy = self.dem.rio.resolution()
        dzdx = (self.dem.shift(x=1) - self.dem.shift(x=-1)) / (2 * dx)
        dzdy = (self.dem.shift(y=1) - self.dem.shift(y=-1)) / (2 * dy)

        slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))[0]

        return slope

    # =====================================================
    # HYDROLOGIC AND TOPOGRAPHIC ANALYSIS
    # =====================================================

    def _calculate_flow_accumulated(self) -> xr.DataArray:
        """
        Calculate flow accumulation from DEM.

        Returns:            Flow accumulation as DataArray
        """
        grid = Grid.from_raster(str(self.dem_path))
        dem = grid.read_raster(str(self.dem_path))

        # Preprocess DEM
        pit_filled_dem = grid.fill_pits(dem)
        flooded_dem = grid.fill_depressions(pit_filled_dem)
        inflated_dem = grid.resolve_flats(flooded_dem)

        # Calculate flow direction and accumulation
        fdir = grid.flowdir(inflated_dem)
        accum = grid.accumulation(fdir)

        accum = xr.DataArray(accum, coords=[self.dem.coords["y"], self.dem.coords["x"]])

        # Convert to area units
        cell_size = self.dem.rio.resolution()[0]
        accum *= cell_size**2

        return accum

    def _calculate_hillshade(
        self, sun_elevation: float = 45, sun_azimuth: float = 315
    ) -> np.ndarray:
        """
        Calculate hillshade from DEM.

        Args:
            sun_elevation: Solar elevation angle in degrees
            sun_azimuth: Solar azimuth angle in degrees

        Returns:
            Hillshade as numpy array
        """
        dem_data = self.dem.values[0]

        sun_elevation_rad = np.radians(sun_elevation)
        sun_azimuth_rad = np.radians(sun_azimuth)

        # Calculate slope and aspect
        gradient_y, gradient_x = np.gradient(dem_data)
        slope = np.arctan(np.sqrt(gradient_x**2 + gradient_y**2))
        aspect = np.arctan2(-gradient_x, gradient_y)

        # Calculate hillshade
        hillshade = (
            255
            * (
                (np.cos(sun_elevation_rad) * np.cos(slope))
                + (
                    np.sin(sun_elevation_rad)
                    * np.sin(slope)
                    * np.cos(aspect - sun_azimuth_rad)
                )
            )
            / 2
            + 128
        )

        return hillshade

    def create_base_plot(self, title: str) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create base plot with hillshade background.

        Args:
            title: Plot title

        Returns:
            Figure and axes objects
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.imshow(self.hillshade, cmap="gray", extent=self.extent)

        plt.title(f"{title}\n")
        plt.axis("off")
        plt.tight_layout()

        return fig, ax

    # =====================================================
    # SHALSTAB STABILITY ANALYSIS METHODS
    # =====================================================

    def _calculate_unconditional_stability(self) -> np.ndarray:
        """
        Calculate unconditional stability conditions.

        Returns:
            Array representing stability conditions:
                * 1: Unconditionally Stable
                * 2: Unconditionally Unstable
                * NaN: Conditionally stable/unstable
        """
        tan_slope = np.tan(self.slope_rad)
        cos_slope_sq = np.cos(self.slope_rad) ** 2

        # Right side for unstable condition
        unstable_threshold = np.tan(self.friction_rad) + (
            self.cohesion / (self.unit_weight * self._soil_thickness * cos_slope_sq)
        )

        # Right side for stable condition
        stable_threshold = (
            (1 - (self.GRAVITY_WATER / self.unit_weight)) * np.tan(self.friction_rad)
        ) + (self.cohesion / (self.unit_weight * self._soil_thickness * cos_slope_sq))

        stability = np.where(
            tan_slope >= unstable_threshold,
            2,  # Unconditionally unstable
            np.where(tan_slope < stable_threshold, 1, np.nan),  # Unconditionally stable
        )

        return stability

    def _plot_unconditional_stability(self, ax: plt.Axes) -> None:
        """Plot unconditional stability classes on given axes."""
        colors = {1: "none", 2: "black"}
        present_values = [v for v in colors.keys() if v in self.unconditional_stability]

        if present_values:
            cmap = ListedColormap([colors[val] for val in present_values])
            ax.imshow(
                self.unconditional_stability,
                cmap=cmap,
                interpolation="nearest",
                extent=self.extent,
                alpha=0.7,
            )  # Add legend
            color_labels = {
                1: "Unconditionally Stable",
                2: "Unconditionally Unstable",
            }

            legend_patches = [
                mpatches.Patch(color=colors[val], label=color_labels[val])
                for val in present_values
            ]
            plt.legend(handles=legend_patches, title="Stability", loc="lower right")

    def calculate_critical_rainfall(self, show_plot: bool = True) -> xr.DataArray:
        """
        Calculate critical rainfall threshold for slope instability using SHALSTAB model.

        This method implements the SHALSTAB equation to determine the rainfall intensity
        required to trigger shallow landslides. The calculation combines topographic,
        hydrologic, and geotechnical factors in a physically-based stability analysis.

        The critical rainfall represents the steady-state precipitation rate that reduces
        the factor of safety to exactly 1.0 (incipient failure condition).

        Parameters
        ----------
        show_plot : bool, default True
            Whether to generate and display a visualization of the results.
            If True, creates a plot with hillshade background, stability classes,
            and critical rainfall values with categorical color scheme.

        Returns
        -------
        xarray.DataArray
            Critical rainfall raster with the following properties:
            - Values: Critical rainfall intensity in mm/day (float64)
            - Coordinates: Same spatial grid as input DEM
            - CRS: Inherits coordinate reference system from DEM
            - NoData: Areas with unconditional stability are set to NaN
            - Units: mm/day (millimeters per day)
            - Range: Typically 0.1 to 1000+ mm/day depending on local conditions

        Raises
        ------
        AttributeError
            If geotechnical parameters have not been initialized (geo_columns not provided
            during class initialization).
        ValueError
            If any of the required geotechnical parameter arrays contain invalid values
            (negative permeability, zero soil thickness, etc.).
        """
        cell_size = self.dem.rio.resolution()[0]
        cos_slope = np.cos(self.slope_rad)
        sin_slope = np.sin(self.slope_rad)
        tan_friction = np.tan(self.friction_rad)

        # Calculate critical rainfall using SHALSTAB equation
        critical = (
            1000
            * (
                (self.permeability)
                * 86400
                * self._soil_thickness
                * cos_slope
                * sin_slope
                * cell_size
                / self.flow_accumulated
            )
            * (
                (
                    self.cohesion
                    / (
                        self.GRAVITY_WATER
                        * self._soil_thickness
                        * cos_slope**2
                        * tan_friction
                    )
                )
                + (
                    (self.unit_weight / self.GRAVITY_WATER)
                    * (1 - (np.tan(self.slope_rad) / tan_friction))
                )
            )
        )

        # Mask unconditional areas
        critical = np.where(self.unconditional_stability > 0, np.nan, critical)
        critical = xr.DataArray(
            critical, coords=[self.dem.coords["y"], self.dem.coords["x"]]
        )

        if show_plot:
            self._plot_critical_rainfall(critical)

        return critical

    def _plot_critical_rainfall(self, critical: xr.DataArray) -> None:
        """Plot critical rainfall results."""
        title = (
            "Physical modeling for landslides\n" "Critical rainfall [mm/day]\nSHALSTAB"
        )

        fig, ax = self.create_base_plot(title)

        # Plot unconditional stability classes
        self._plot_unconditional_stability(ax)

        # Plot critical rainfall with custom colormap
        valid_critical = critical.where(critical > 0)
        # Use percentiles on masked DataArray directly as in original
        percentiles = [1, 20, 40, 60, 80, 99]
        colors = np.array([np.nanpercentile(valid_critical, p) for p in percentiles])

        cmap = ListedColormap(self.CRITICAL_COLORS)
        norm = BoundaryNorm(colors, cmap.N)

        plt.imshow(
            critical,
            cmap=cmap,
            interpolation="nearest",
            extent=self.extent,
            alpha=0.7,
            norm=norm,
        )
        plt.colorbar()

    def calculate_log_qt(self) -> xr.DataArray:
        """
        Calculate logarithm of the ratio between steady-state recharge (q) and soil transmissivity (T).

        Returns:
            Log(q/T) raster
        """
        cell_size = self.dem.rio.resolution()[0]
        cos_slope = np.cos(self.slope_rad)
        sin_slope = np.sin(self.slope_rad)
        tan_friction = np.tan(self.friction_rad)

        log_qt = np.log(
            (sin_slope * cell_size / self.flow_accumulated)
            * (
                (
                    self.cohesion
                    / (
                        self.GRAVITY_WATER
                        * self._soil_thickness
                        * cos_slope**2
                        * tan_friction
                    )
                )
                + (
                    (
                        self.unit_weight / self.GRAVITY_WATER
                    )  # original geohazards uses inverse ratio
                    * (1 - (np.tan(self.slope_rad) / tan_friction))
                )
            )
        )

        # Mask unconditional areas
        log_qt = np.where(self.unconditional_stability > 0, np.nan, log_qt)
        log_qt = xr.DataArray(
            log_qt, coords=[self.dem.coords["y"], self.dem.coords["x"]]
        )

        log_qt.rio.write_nodata(self.NODATA_VALUE, inplace=True)
        log_qt.rio.write_crs(self.dem.rio.crs, inplace=True)

        self._plot_log_qt(log_qt)
        return log_qt

    def _plot_log_qt(self, log_qt: xr.DataArray) -> None:
        """Plot log(q/T) results."""
        title = "Physical modeling for landslides\nlog (q/T)\nSHALSTAB"
        fig, ax = self.create_base_plot(title)

        # Plot unconditional stability
        self._plot_unconditional_stability(ax)

        # Plot log(q/T) with custom colormap
        percentiles = [1, 20, 40, 60, 80, 99]
        colors = np.array([np.nanpercentile(log_qt, p) for p in percentiles])

        cmap = ListedColormap(self.CRITICAL_COLORS)
        norm = BoundaryNorm(colors, cmap.N)

        plt.imshow(
            log_qt,
            cmap=cmap,
            interpolation="nearest",
            extent=self.extent,
            alpha=0.7,
            norm=norm,
        )

        cbar = plt.colorbar()
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    def calculate_stability(
        self, rainfall_mm_day: float
    ) -> Tuple[xr.DataArray, plt.Figure]:
        """
        Calculate slope stability classification for a specified rainfall intensity.

        This method applies the SHALSTAB model to classify each pixel into stability
        categories based on infinite slope stability analysis under steady-state
        hydrologic conditions. The analysis compares the stabilizing forces
        (cohesion and friction) against destabilizing forces (gravity and pore pressure).

        Parameters
        ----------
        rainfall_mm_day : float
            Rainfall intensity for stability analysis in millimeters per day.
            Must be positive. Typical values range from 1 to 500 mm/day.
            Higher values represent more intense precipitation events.

        Returns
        -------
        Tuple[xarray.DataArray, matplotlib.figure.Figure]
            A tuple containing:

            stability_raster : xarray.DataArray
                Stability classification grid with integer codes:
                - 1: Unconditionally Stable (stable regardless of rainfall)
                - 2: Unconditionally Unstable (unstable regardless of rainfall)
                - 3: Unstable (unstable at specified rainfall intensity)
                - 4: Stable (stable at specified rainfall intensity)

                Properties:
                - dtype: int32
                - Coordinates: Same as input DEM
                - CRS: Inherits from DEM
                - Attributes: Contains statistical report as 'Reporte' attribute

            figure : matplotlib.figure.Figure
                Visualization of stability results with:
                - Hillshade background for terrain context
                - Categorical color scheme for stability classes
                - Legend with class descriptions
                - Professional formatting and labels

        Raises
        ------
        ValueError
            If rainfall_mm_day is not positive or if geotechnical parameters
            are not initialized.
        AttributeError
            If required attributes (soil properties, topographic derivatives)
            are missing due to incomplete initialization.
        """
        cell_size = self.dem.rio.resolution()[0]
        cos_slope = np.cos(self.slope_rad)
        sin_slope = np.sin(self.slope_rad)
        tan_friction = np.tan(self.friction_rad)

        # Left side of SHALSTAB equation
        left_side = self.flow_accumulated / cell_size

        # Right side of SHALSTAB equation
        right_side = (
            (
                (self.permeability)
                * 86400
                * (self._soil_thickness * cos_slope)
                * sin_slope
            )
            / (rainfall_mm_day / 1000)
        ) * (
            (
                self.GRAVITY_WATER / self.unit_weight
            )  # original geohazards uses inverse ratio
            * (1 - np.tan(self.slope_rad) / tan_friction)
            + (
                self.cohesion
                / (
                    self.GRAVITY_WATER
                    * self._soil_thickness
                    * cos_slope**2
                    * tan_friction
                )
            )
        )

        # Calculate stability categories
        stability = np.where(
            self.unconditional_stability > 0,
            self.unconditional_stability,  # Keep unconditional values
            np.where(left_side > right_side, 3, 4),  # Unstable  # Stable
        )

        stability = xr.DataArray(
            stability, coords=[self.dem.coords["y"], self.dem.coords["x"]]
        )

        # Generate report
        report = self._generate_stability_report(stability)
        stability.attrs["Reporte"] = report

        # Create plot
        fig = self._plot_stability(stability, rainfall_mm_day)

        return stability, fig

    def _generate_stability_report(self, stability: xr.DataArray) -> str:
        """Generate stability statistics report."""
        valid_cells = stability.count().item()
        category_counts = stability.to_series().value_counts()

        report = ""
        for category, label in self.STABILITY_CATEGORIES.items():
            count = category_counts.get(category, 0)
            percentage = count * 100 / valid_cells if valid_cells > 0 else 0
            report += f"{label}: {percentage:.2f}%\n"

        return report

    def _plot_stability(self, stability: xr.DataArray, rainfall: float) -> plt.Figure:
        """Plot stability results."""
        title = (
            f"Physical modeling for landslides\n"
            f"Stability for rainfall of {rainfall} mm/day\nSHALSTAB"
        )

        fig, ax = self.create_base_plot(title)

        # Get present stability values and their colors
        present_values = [
            val for val in self.STABILITY_COLORS.keys() if val in stability.values
        ]
        present_colors = [self.STABILITY_COLORS[val] for val in present_values]
        present_labels = [self.STABILITY_CATEGORIES[val] for val in present_values]

        if present_values:
            cmap = ListedColormap(present_colors)
            plt.imshow(
                stability,
                cmap=cmap,
                interpolation="nearest",
                extent=self.extent,
                alpha=0.7,
            )  # Create legend
            legend_patches = [
                mpatches.Patch(color=color, label=label)
                for color, label in zip(present_colors, present_labels)
            ]

            plt.legend(
                handles=legend_patches,
                title="Stability",
                loc="lower center",
                bbox_to_anchor=(0.5, -0.2),
                ncol=2,
            )

        return fig

    def calculate_failure_probability(self) -> xr.DataArray:
        """
        Calculate relative failure probability based on topographic and hydrologic factors.

        This method derives a relative probability of landslide failure (0-100%) from
        the log(q/T) ratio calculated by the SHALSTAB model. The probability represents
        the relative likelihood of failure across the landscape, with higher values
        indicating areas more susceptible to shallow landsliding.

        The calculation normalizes the inverted log(q/T) values to create a probability
        surface where areas requiring lower critical rainfall (more unstable) receive
        higher probability scores.

        Returns
        -------
        xarray.DataArray
            Relative failure probability raster with the following properties:

            Values and Range:
            - Data type: float64
            - Range: 0.0 to 100.0 (percentage)
            - Units: Percent probability (%)
            - Higher values = greater failure likelihood
            - Lower values = lower failure likelihood

            Spatial Properties:
            - Coordinates: Same spatial grid as input DEM
            - CRS: Inherits coordinate reference system from DEM
            - Resolution: Same as input DEM
            - NoData: Areas with unconditional stability set to NaN            Attributes:
            - 'Reporte': Statistical summary with percentage distribution
              across probability ranges (0-20%, 20-40%, 40-60%, 60-80%, 80-100%)
            - CRS information and NoData value metadata

        Raises
        ------
        AttributeError
            If geotechnical parameters are not initialized or calculate_log_qt()
            cannot be executed due to missing required data.
        RuntimeError
            If log(q/T) calculation fails due to numerical issues (infinite values,
            all NaN results, etc.).
        """
        log_qt = self.calculate_log_qt()
        log_qt_inverted = -log_qt

        # Normalize to 0-100% probability
        prob = (
            (log_qt_inverted - log_qt_inverted.min())
            / (log_qt_inverted.max() - log_qt_inverted.min())
            * 100
        )

        # Mask unconditional areas
        prob = np.where(self.unconditional_stability > 0, np.nan, prob)
        prob = xr.DataArray(prob, coords=[self.dem.coords["y"], self.dem.coords["x"]])

        prob.rio.write_crs(self.dem.rio.crs, inplace=True)
        prob.rio.write_nodata(self.NODATA_VALUE, inplace=True)

        # Generate report
        report = self._generate_probability_report(prob)
        prob.attrs["Reporte"] = report

        # Create plot
        self._plot_probability(prob)

        return prob

    def _generate_probability_report(self, prob: xr.DataArray) -> str:
        """Generate probability statistics report."""
        total_cells = prob.count().item()
        report = ""

        ranges = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
        range_labels = ["0-20", "20-40", "40-60", "60-80", "80-100"]

        for (low, high), label in zip(ranges, range_labels):
            range_cells = np.count_nonzero((prob >= low) & (prob < high))
            percentage = range_cells * 100 / total_cells if total_cells > 0 else 0
            report += f"{label}%: {percentage:.2f}%\n"

        return report

    def _plot_probability(self, prob: xr.DataArray) -> None:
        """Plot failure probability results."""
        title = "Physical modeling for landslides\n" "Failure probability\nSHALSTAB"

        fig, ax = self.create_base_plot(title)

        # Plot unconditional stability
        self._plot_unconditional_stability(ax)  # Plot probability with custom colormap
        colors = np.array([0, 20, 40, 60, 80, 100])
        cmap = ListedColormap(self.PROBABILITY_COLORS)
        norm = BoundaryNorm(colors, cmap.N)

        plt.imshow(
            prob,
            cmap=cmap,
            interpolation="nearest",
            extent=self.extent,
            alpha=0.7,
            norm=norm,
        )
        plt.colorbar()
