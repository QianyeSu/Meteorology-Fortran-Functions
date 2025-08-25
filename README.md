# Fortran Code Function Analysis

## 1. KolmSmir2.f [📝 src](fortran/KolmSmir2.f)
**Function**: Kolmogorov-Smirnov two-sample test implementation

**Detailed Description**:
- Used to test differences between two sample distribution functions
- Contains two subroutines:
  - KOLM2: Main test routine that performs statistical testing on input samples X and Y
  - SMIRN: Calculates probability values corresponding to statistics
- Input parameters:
  - X: Sample 1 data vector (length N)
  - Y: Sample 2 data vector (length M)
  - N,M: Size of the two samples
  - IFLAG: Control flag
- Output parameters:
  - Z: Test statistic
  - PROB: Corresponding probability value
  - D: Maximum difference between the two distribution functions
- Key features:
  - Automatically sorts input data in ascending order
  - Calculates maximum difference between empirical distribution functions
  - Suitable for large sample testing (N,M>=100)
  - Based on classical Kolmogorov-Smirnov statistical theory

## 2. VSPCORR.f [📝 src](fortran/VSPCORR.f)
**Function**: Spearman rank correlation coefficient calculation (missing value handling version)

**Detailed Description**:
- Contains four subroutines:
  - SPCORRZ: Main entry program that handles data with missing values
  - SPCORR: Core program for calculating Spearman rank correlation coefficient
  - SPRANK: Data ranking routine
  - SPSORT: Binary sorting algorithm
- Input parameters:
  - X,Y: Input data vectors (length N)
  - N: Number of samples
  - HASXMSG,HASYMSG: Missing value flags
  - XMSG,YMSG: Missing value codes
  - IWRITE: Output control parameter
- Output parameters:
  - SPC: Spearman rank correlation coefficient (between -1 and 1)
  - NMSG: Number of missing values
- Key features:
  - Automatically handles and excludes missing values
  - Uses efficient binary sorting algorithm
  - Supports average rank handling for tied values
  - Suitable for non-parametric correlation analysis
  - Based on classical Spearman rank correlation theory

## 3. aam.f [📝 src](fortran/aam.f)
**Function**: Atmospheric Angular Momentum calculation

**Detailed Description**:
- Contains three subroutines:
  - DAAMOM1: Calculate atmospheric angular momentum (using 1D pressure thickness array)
  - DAAMOM3: Calculate atmospheric angular momentum (using 3D pressure thickness array)
  - DCOSWEIGHT: Calculate latitude weights
- Input parameters:
  - U: Zonal wind component (m/s), 3D array (longitude×latitude×levels)
  - DP: Pressure thickness (Pa), can be 1D or 3D array
  - LAT: Latitude array
  - WGT: Latitude weight array
  - MLON,NLAT,KLVL: Number of longitude, latitude, and vertical levels
  - UMSG: Wind speed missing value code
- Output parameters:
  - AAM: Atmospheric relative angular momentum (kg·m²/s)
- Key features:
  - Based on geophysical angular momentum conservation principle
  - Supports different dimensional pressure thickness inputs
  - Automatically handles missing values
  - Considers Earth radius and gravitational acceleration
  - Uses cosine latitude weighting
  - Suitable for climate and atmospheric circulation analysis

## 4. accumrun.f [📝 src](fortran/accumrun.f)
**Function**: Running window accumulation (Running Accumulation)

**Detailed Description**:
- Contains one subroutine:
  - DACUMRUN: Calculate sliding accumulation with specified window size
- Input parameters:
  - P: Input time series data array (length NP)
  - NP: Total number of data points
  - PMSG: Missing value code
  - NRUN: Sliding window size
  - IOPT: Missing value handling option (0 or 1)
- Output parameters:
  - PAC: Accumulation array
  - NMSG: Total number of missing values encountered
- Key features:
  - First (NRUN-1) values are set to missing
  - Starts calculating sliding accumulation from the NRUN-th point
  - Supports two missing value handling strategies:
    - IOPT=0: Terminates current window calculation when encountering missing values, sets result to missing
    - IOPT=1: Counts missing values but continues calculating remaining valid values
  - Suitable for precipitation accumulation, temperature accumulation, and other time series analysis

## 5. areaAve.f [📝 src](fortran/areaAve.f)
**Function**: Calculate weighted area average of 2D data (using separated weights)

**Detailed Description**:
- Contains one subroutine:
  - DWGTAREAAVE: Calculate area average using separated X and Y direction weights
- Input parameters:
  - T: 2D data array (MX×NY)
  - WGTX: X-direction weight array (length MX, e.g., longitude weights)
  - WGTY: Y-direction weight array (length NY, e.g., latitude weights)
  - MX,NY: Array dimensions (longitude×latitude)
  - TMSG: Missing value code
  - IFLAG: Missing value handling flag (0=ignore missing values, 1=return missing value when encountered)
- Output parameters:
  - AVE: Weighted area average
- Key features:
  - Supports area weighting for longitude-latitude grids
  - Weights are the product of WGTX(i)×WGTY(j)
  - Flexible missing value handling strategy
  - Suitable for spatial averaging calculations of climate data

## 6. areaAve2.f [📝 src](fortran/areaAve2.f)
**Function**: Calculate weighted area average of 2D data (using 2D weight matrix)

**Detailed Description**:
- Contains one subroutine:
  - DWGTAREAAVE2: Calculate area average using 2D weight matrix
- Input parameters:
  - T: 2D data array (MX×NY)
  - WGT: 2D weight matrix (MX×NY)
  - MX,NY: Array dimensions
  - TMSG: Missing value code
  - IFLAG: Missing value handling flag
- Output parameters:
  - AVE: Weighted area average
- Key features:
  - Supports arbitrarily complex weight distributions
  - Each grid point can have independent weight values
  - Suitable for irregular grids or special weighting requirements
  - Flexible missing value handling

## 7. areaRmse.f [📝 src](fortran/areaRmse.f)
**Function**: Calculate weighted root mean square error between two 2D data arrays (using separated weights)

**Detailed Description**:
- Contains one subroutine:
  - DWGTAREARMSE: Calculate RMSE using separated X and Y direction weights
- Input parameters:
  - T,Q: Two 2D data arrays (MX×NY) for comparison
  - WGTX: X-direction weight array (length MX)
  - WGTY: Y-direction weight array (length NY)
  - MX,NY: Array dimensions
  - TMSG,QMSG: Missing value codes for T and Q
  - IFLAG: Missing value handling flag
- Output parameters:
  - RMSE: Weighted root mean square error
- Key features:
  - Calculation formula: $\text{RMSE} = \sqrt{\frac{\sum \text{wgt}(T-Q)^2}{\sum \text{wgt}}}$
  - Supports error analysis between two datasets
  - Weights are the product of $\text{WGTX}(i) \times \text{WGTY}(j)$
  - Suitable for model validation and data comparison

## 8. areaRmse2.f [📝 src](fortran/areaRmse2.f)
**Function**: Calculate weighted root mean square error between two 2D data arrays (using 2D weight matrix)

**Detailed Description**:
- Contains one subroutine:
  - DWGTAREARMSE2: Calculate RMSE using 2D weight matrix
- Input parameters:
  - T,Q: Two 2D data arrays (MX×NY)
  - WGT: 2D weight matrix (MX×NY)
  - MX,NY: Array dimensions
  - TMSG,QMSG: Missing value codes
  - IFLAG: Missing value handling flag
- Output parameters:
  - RMSE: Weighted root mean square error
- Key features:
  - Supports arbitrarily complex weight distributions
  - Each grid point can have independent weights
  - Suitable for error analysis on irregular grids
  - Flexible missing value handling

## 9. areaSum2.f [📝 src](fortran/areaSum2.f)
**Function**: Calculate weighted area sum of 2D data (using 2D weight matrix)

**Detailed Description**:
- Contains one subroutine:
  - DWGTAREASUM2: Calculate area weighted sum using 2D weight matrix
- Input parameters:
  - T: 2D data array (MX×NY)
  - WGT: 2D weight matrix (MX×NY)
  - MX,NY: Array dimensions
  - TMSG: Missing value code
  - IFLAG: Missing value handling flag
- Output parameters:
  - SUM: Weighted area sum
- Key features:
  - Calculation formula: $\text{SUM} = \sum_{i=1}^{MX} \sum_{j=1}^{NY} T(i,j) \times \text{WGT}(i,j)$
  - Supports arbitrarily complex weight distributions
  - Each grid point can have independent weight values
  - Flexible missing value handling strategies:
    - IFLAG=0: Ignore missing values, continue calculating sum of valid values
    - IFLAG=1: Return missing value immediately when encountering missing values
  - Suitable for area integration calculations on irregular grids
  - Flexible handling of missing value situations

## 10. area_conremap.f [📝 src](fortran/area_conremap.f)
**Function**: Conservative remapping for grid point interpolation

**Detailed Description**:
- Contains four subroutines:
  - CREMAPBIN: Main remapping program using grid box binning method
  - BINNING_GET_GLOBAL_LATS_WGTS: Get global Gaussian latitudes and weights
  - BINNING_MAP_LATS: Map regional latitudes to global latitudes
  - BINNING_MAP_EDGES: Calculate grid box boundaries
- Input parameters:
  - XX: Input 3D data field (longitude×latitude×levels)
  - CLAT,CLON: Input grid latitudes and longitudes
  - CLATO,CLONO: Output grid latitudes and longitudes
  - BIN_FACTOR: Binning area expansion/contraction factor
  - NLAT,NLATO: Number of global Gaussian latitudes
- Output parameters:
  - YY: Output field after horizontal interpolation
- Key features:
  - Supports both regular and Gaussian grids
  - Conservative interpolation algorithm maintaining mass conservation
  - Supports vertically multi-level data
  - Automatically handles missing values
  - Suitable for climate model data regridding

## 11. bfilter1.f [📝 src](fortran/bfilter1.f)
**Function**: Butterworth bandpass filter

**Detailed Description**:
- Contains two subroutines:
  - BUTTFILT: Filter entry program
  - BFILTER: Core filter implementation
- Input parameters:
  - XR: Input time series (length N)
  - M: Butterworth filter order
  - FCA,FCB: Low and high frequency cutoff frequencies
  - DT: Time interval
  - MZER: Flag for mean removal
- Output parameters:
  - YR: Filtered time series
  - ER: Envelope function
  - IER: Error code
- Key features:
  - Zero-phase Butterworth bandpass filtering
  - Cascaded first-order filter implementation
  - Forward and backward filtering ensures zero phase
  - Complex time series processing
  - Stable filtering algorithm
  - Suitable for seismology and signal processing

## 12. binavg.f [📝 src](fortran/binavg.f)
**Function**: Grid binning average for irregularly distributed data

**Detailed Description**:
- Contains one subroutine:
  - BINDATAAVG: Grid averaging of irregular geographic data
- Input parameters:
  - ZLON,ZLAT: Longitude and latitude coordinates of data points (NZ points)
  - Z: Data value array (length NZ)
  - GLON,GLAT: Target grid longitude and latitude (MLON×NLAT)
  - ZMSG: Missing value identifier
  - IOPT: Option parameter
- Output parameters:
  - GBINKNT: Grid binning results (1st layer: averages, 2nd layer: counts)
  - IER: Error code
- Key features:
  - Validates grid equal spacing
  - Assigns data points to corresponding grid cells
  - Calculates average value within each grid cell
  - Suitable for geographic data regridding

## 13. bindata.f [📝 src](fortran/bindata.f)
**Function**: Grid binning summation for irregularly distributed data

**Detailed Description**:
- Contains one subroutine:
  - BINDATASUM3: Grid summation of irregular geographic data
- Input parameters:
  - ZLON,ZLAT: Data point coordinates
  - Z: Data value array
  - GLON,GLAT: Target grid coordinates
  - ZMSG: Missing value identifier
- Output parameters:
  - GBIN: Grid binning summation results
  - GKNT: Data point count for each grid cell
- Key features:
  - Performs geographic data grid binning summation
  - Provides data redistribution statistics
  - Preprocessing for statistical analysis
  - No averaging calculation, only summation

## 14. calc_uh.f90 [📝 src](fortran/calc_uh.f90)
**Function**: Calculate Updraft Helicity

**Detailed Description**:
- Contains one subroutine:
  - DCALCUH: Calculate UH diagnostic quantity
- Input parameters:
  - US,VS: u and v wind components (nx×ny×nz)
  - W: Vertical velocity (nx×ny×nzp1)
  - ZP: Height array
  - MAPFCT: Map projection factor
  - UHMNHGT,UHMXHGT: Calculation height range
- Output parameters:
  - UH: Updraft helicity field (nx×ny)
- Key features:
  - Based on Kain et al. 2008 formulation
  - Integrates vertical vorticity × vertical velocity over specified height range
  - Used for identifying supercell storms and tornadoes
  - Uses OpenMP parallel computation
  - Important meteorological diagnostic quantity

## 15. calcorc_dp.f [📝 src](fortran/calcorc_dp.f)
**Function**: Calculate Pearson correlation coefficient and covariance

**Detailed Description**:
- Contains two subroutines:
  - DCALCORC: Calculate correlation coefficient
  - DCALCOVC: Calculate covariance
- Input parameters:
  - X,Y: Two data sequences (length NXY)
  - XAVE,XSTD: Mean and standard deviation of X sequence
  - XMSG,YMSG: Missing value identifiers
- Output parameters:
  - CORC: Pearson correlation coefficient
  - COVC: Covariance
  - IER: Error code
- Key features:
  - Capable of handling missing values
  - Automatic statistical quantity calculation
  - Standardized and unstandardized calculations
  - Fundamental tool for statistical analysis

## 16. cancor.f [📝 src](fortran/cancor.f)
**Function**: Canonical Correlation Analysis

**Detailed Description**:
- Contains multiple subroutines:
  - DCANCORXY: Main analysis program
  - DCANORX,DMINVX,DNROOTX,DEIGENX: Matrix operation subroutines
- Input parameters:
  - X: Independent variable matrix (nobs×mx)
  - Y: Dependent variable matrix (nobs×my)
  - NOBS: Number of observations
- Output parameters:
  - CANR: Canonical correlation coefficients
  - EVAL: Eigenvalues
  - CHISQ: Chi-square statistics
  - COEFXR,COEFYL: Canonical coefficients
- Key features:
  - Identifies linear relationship patterns between two variable sets
  - Contains complete matrix operation suite
  - Advanced multivariate statistical analysis technique
  - Used for dimensionality reduction and pattern recognition

## 17. cdf_dp.f [📝 src](fortran/cdf_dp.f)
**Function**: Cumulative distribution function calculation for multiple probability distributions

**Detailed Description**:
- Contains multiple subroutines:
  - DCDFBINP,DCDFBINS,DCDFBINXN: Binomial distribution CDF
  - DCDFGAMP,DCDFGAMX: Gamma distribution CDF
  - DCDFNORP,DCDFNORX: Normal distribution CDF
  - DCDFCHIP: Chi-square distribution CDF
- Input parameters: Distribution parameters (probability, degrees of freedom, shape parameters, etc.)
- Output parameters: Corresponding CDF values or inverse function values
- Key features:
  - Supports multiple probability distributions
  - Provides forward and inverse calculations
  - Comprehensive error handling
  - Core tool for probability statistics and hypothesis testing

## 18. cdft_dp.f [📝 src](fortran/cdft_dp.f)
**Function**: Cumulative distribution function calculation for t-distribution

**Detailed Description**:
- Contains three subroutines:
  - DCDFCDFTT,DCDFCDFTP: t-distribution CDF calculation
  - DCDFTERR: Error handling
- Input parameters:
  - P: Probability value array
  - T: t-statistic array
  - DF: Degrees of freedom array
- Output parameters:
  - Output t-values or probability values depending on calculation direction
  - IER: Error code
- Key features:
  - Specialized for t-distribution calculations
  - Provides bidirectional calculation functionality
  - Comprehensive error handling mechanism
  - Fundamental tool for t-tests and confidence interval calculations

## 19. cfft_driver.f [📝 src](fortran/cfft_driver.f)
**Function**: Complex FFT driver program

**Detailed Description**:
- Contains four subroutines:
  - CFFTFDRIVER: Complex forward FFT driver
  - CFFTBDRIVER: Complex inverse FFT driver
  - FRQCFFT: Generate FFT frequencies
  - CFFTFFRQREORDER: Reorder frequencies
- Input parameters: Complex real and imaginary part arrays, work arrays
- Output parameters: FFT coefficients, frequency arrays
- Key features:
  - Provides high-level complex FFT interface for NCL
  - Includes forward and inverse transform functionality
  - Supports frequency generation and reordering
  - Suitable for time series spectral analysis and signal processing

## 20. cfftb.f [📝 src](fortran/cfftb.f)
**Function**: Complex FFT inverse transform

**Detailed Description**:
- Contains one subroutine:
  - CFFTB: Complex discrete Fourier inverse transform
- Input parameters: Complex sequence C, work array WSAVE
- Output parameters: Inverse transformed complex sequence
- Key features:
  - Performs Fourier synthesis
  - Reconstructs complex periodic sequence from Fourier coefficients
  - Frequency domain to time domain conversion
  - Used for signal reconstruction and post-filtering recovery

## 21. cfftb1.f [📝 src](fortran/cfftb1.f)
**Function**: Complex FFT inverse transform core algorithm

**Detailed Description**:
- Contains one subroutine:
  - CFFTB1: Complex FFT inverse transform core implementation
- Input parameters: Sequence length N, complex array C, work arrays
- Output parameters: Transform results
- Key features:
  - Low-level inverse transform implementation
  - Handles FFT calculations with different factor decompositions
  - Supports basic factors like 2, 3, 4, 5, etc.
  - Internal core computation of FFT algorithm

## 22. cfftf.f [📝 src](fortran/cfftf.f)
**Function**: Complex FFT forward transform

**Detailed Description**:
- Contains one subroutine:
  - CFFTF: Complex discrete Fourier forward transform
- Input parameters: Complex sequence C, work array WSAVE
- Output parameters: Fourier coefficients
- Key features:
  - Performs Fourier analysis
  - Calculates Fourier coefficients of complex periodic sequences
  - Time domain to frequency domain conversion
  - Used for spectral analysis and signal decomposition

## 23. cfftf1.f [📝 src](fortran/cfftf1.f)
**Function**: Complex FFT forward transform core algorithm

**Detailed Description**:
- Contains one subroutine:
  - CFFTF1: Complex FFT forward transform core implementation
- Input parameters: Sequence length N, complex array C, work arrays
- Output parameters: Transform results
- Key features:
  - Low-level forward transform implementation
  - Uses mixed-radix algorithm
  - Handles sequences of different lengths
  - Efficient forward transform processing

## 24. cffti.f [📝 src](fortran/cffti.f)
**Function**: Complex FFT initialization

**Detailed Description**:
- Contains one subroutine:
  - CFFTI: Complex FFT initialization program
- Input parameters: Sequence length N
- Output parameters: Work array WSAVE
- Key features:
  - Initializes work arrays required for FFT computation
  - Calculates and stores prime factor decomposition
  - Generates trigonometric function tables
  - FFT preprocessing and work array preparation

## 25. cffti1.f [📝 src](fortran/cffti1.f)
**Function**: Complex FFT initialization core algorithm

**Detailed Description**:
- Contains one subroutine:
  - CFFTI1: Complex FFT initialization core implementation
- Input parameters: Sequence length N
- Output parameters: Work array WA, factor array IFAC
- Key features:
  - Performs actual FFT initialization work
  - Prime factor decomposition
  - Trigonometric function calculation and storage
  - Optimized precomputation for repeated transforms

## 26. chiin_dp.f [📝 src](fortran/chiin_dp.f)
**Function**: Chi-square distribution inverse function calculation

**Detailed Description**:
- Contains one subroutine:
  - CHISUB: Chi-square distribution inverse function interface
- Input parameters: Cumulative probability P, degrees of freedom DF
- Output parameters: Chi-square value CHI
- Key features:
  - Calculates inverse function of chi-square distribution
  - Given probability and degrees of freedom, finds chi-square statistic
  - Suitable for statistical hypothesis testing
  - Confidence interval calculation and probability analysis

## 27. cntrFiniteDiff.f [📝 src](fortran/cntrFiniteDiff.f)
**Function**: Central finite difference calculation

**Detailed Description**:
- Contains one subroutine:
  - DCFINDIF: Central finite difference calculation
- Input parameters: Variable array Q, coordinate array R, boundary condition parameters
- Output parameters: Derivative array DQDR
- Key features:
  - Uses central finite difference method
  - Calculates numerical derivatives
  - Supports cyclic and non-cyclic boundary conditions
  - Used for numerical differentiation and gradient calculations

## 28. conservLinInt_2d_ncl.f [📝 src](fortran/conservLinInt_2d_ncl.f)
**Function**: Two-dimensional conservative linear interpolation

**Detailed Description**:
- Contains one subroutine:
  - AREALININT2DA: Two-dimensional area-weighted interpolation
- Input parameters: Input grid coordinates, data, weights
- Output parameters: Interpolated grid data
- Key features:
  - Implements two-dimensional conservative linear interpolation
  - Area-weighted method
  - Maintains integral conservation of data
  - Suitable for grid remapping and climate data processing

## 29. covcorm_driver.f [📝 src](fortran/covcorm_driver.f)
**Function**: Covariance correlation matrix calculation driver

**Detailed Description**:
- Contains two subroutines:
  - DCOVCORMSSM: Symmetric storage mode covariance matrix
  - DCOVCORM: Two-dimensional covariance matrix
- Input parameters: Data matrix, missing value identifiers, options
- Output parameters: Covariance or correlation matrix
- Key features:
  - Calculates covariance matrix or correlation matrix
  - Supports symmetric storage and full matrix storage
  - Multivariate statistical analysis tool
  - Used for principal component analysis and correlation analysis

## 30. covcorm_xy_matrix_ncl.f [📝 src](fortran/covcorm_xy_matrix_ncl.f)
**Function**: Cross-covariance matrix calculation

**Detailed Description**:
- Contains one subroutine:
  - DCOVARXY: Bivariate covariance calculation
- Input parameters: Variable arrays X, Y, lag parameters
- Output parameters: Cross-covariance/correlation matrix
- Key features:
  - Calculates cross-covariance between two multidimensional variables
  - Supports time lag analysis
  - Cross-correlation analysis tool
  - Suitable for lag correlation studies and time series analysis

## 31. ctwrap.f [📝 src](fortran/ctwrap.f)
**Function**: Three-dimensional grid visualization system

**Detailed Description**:
- Contains multiple subroutines: CTDRIVER, DFCLRS, DRWMCL, DCFOCB, MKMESH, DSJE02, etc.
- Input parameters: Triangular grid data, plotting parameters, work arrays
- Output parameters: Various visualization graphics and triangular grids
- Key features:
  - Complex three-dimensional grid visualization system
  - Generates triangular grids from SEAM meshes
  - Supports multiple visualization methods
  - Used for 3D visualization of Earth surface grid data, contour plotting, and filled contour maps

## 32. cz2ccm_dp.f [📝 src](fortran/cz2ccm_dp.f)
**Function**: CCM2 hybrid coordinate system geophysical height calculation

**Detailed Description**:
- Contains two subroutines: DCZ2CCM, DCZ2
- Input parameters: Surface pressure PS, surface geopotential PHIS, virtual temperature TV, hybrid coordinate coefficients
- Output parameters: Mid-layer geopotential height Z2
- Key features:
  - Based on CCM2 hybrid coordinate system
  - Calculates geophysical height fields
  - Used for vertical coordinate conversion and meteorological data processing in atmospheric models

## 33. d1mach_alpha.f [📝 src](fortran/d1mach_alpha.f)
**Function**: Alpha platform double precision floating-point constants

**Detailed Description**:
- Contains one subroutine: D1MACH
- Input parameters: Integer parameter I (1-5)
- Output parameters: Alpha architecture machine-dependent double precision floating-point constants
- Key features:
  - DEC Alpha processor optimized floating-point constants
  - IEEE 754 standard compatible
  - High-precision numerical computation support

## 34. d1mach_linux.f [📝 src](fortran/d1mach_linux.f)
**Function**: Linux platform double precision floating-point constants

**Detailed Description**:
- Contains one subroutine: D1MACH
- Input parameters: Integer parameter I (1-5)
- Output parameters: Linux system machine-dependent double precision floating-point constants
- Key features:
  - x86/x64 architecture optimization
  - GNU compiler compatibility
  - Cross-Linux distribution support

## 35. d1mach_sun.f [📝 src](fortran/d1mach_sun.f)
**Function**: Sun Solaris platform double precision floating-point constants

**Detailed Description**:
- Contains one subroutine: D1MACH
- Input parameters: Integer parameter I (1-5)
- Output parameters: Sun SPARC architecture machine-dependent double precision floating-point constants
- Key features:
  - SPARC processor optimization
  - Sun compiler compatibility
  - Solaris system integration

## 36. dbetai.f [📝 src](fortran/dbetai.f)
**Function**: Incomplete Beta function calculation

**Detailed Description**:
- Contains main programs: BETAINC, DBETAISLATEC and multiple SLATEC auxiliary functions
- Input parameters: Integration upper limit x, beta distribution parameters a,b
- Output parameters: Incomplete Beta function value result
- Key features:
  - Calculates incomplete Beta function
  - Contains complete SLATEC mathematical library support
  - Used for statistical analysis, probability calculations, hypothesis testing

## 37. dcolconv.f [📝 src](fortran/dcolconv.f)
**Function**: Color space conversion library

**Detailed Description**:
- Contains multiple subroutines: DHLSRGB, DHSVRGB, DRGBHLS, DRGBHSV, DRGBYIQ, DYIQRGB, etc.
- Input parameters: Color values in various color spaces
- Output parameters: Converted color values
- Key features:
  - Complete color space conversion library
  - Supports mutual conversion between RGB, HLS, HSV, YIQ
  - Used for graphics visualization, color processing, image processing

## 38. det_code42.f [📝 src](fortran/det_code42.f)
**Function**: Matrix determinant calculation

**Detailed Description**:
- Contains three subroutines: DETCALC, DTRM, ELGS
- Input parameters: n×n matrix a
- Output parameters: Matrix determinant value det
- Key features:
  - Uses partial pivoting Gaussian elimination
  - Calculates matrix determinant
  - Used for linear algebra calculations, numerical analysis, solving linear equation systems

## 39. dewpt_trh.f [📝 src](fortran/dewpt_trh.f)
**Function**: Dew point temperature calculation

**Detailed Description**:
- Contains one subroutine: DDEWTEMP
- Input parameters: Temperature TK[K], relative humidity RH[%]
- Output parameters: Dew point temperature TDK[K]
- Key features:
  - Calculates dew point temperature from temperature and relative humidity
  - Used for meteorological analysis, humidity calculations, atmospheric physics research

## 40. dgeevx_interface.f [📝 src](fortran/dgeevx_interface.f)
**Function**: LAPACK eigenvalue calculation interface

**Detailed Description**:
- Contains one subroutine: DGEEVXINT
- Input parameters: Matrix A, control parameters, work arrays
- Output parameters: Eigenvalues WR,WI and eigenvectors EVLR
- Key features:
  - Interface to LAPACK DGEEVX function
  - Calculates matrix eigenvalues and eigenvectors
  - Used for numerical linear algebra, principal component analysis, dynamical system analysis

## 41. dimavgwgt.f [📝 src](fortran/dimavgwgt.f)
**Function**: One-dimensional weighted average calculation

**Detailed Description**:
- Contains one subroutine: DIMAVGWGT
- Input parameters: Array x, weights wgt, control parameter iopt
- Output parameters: Weighted average xavg
- Key features:
  - Calculates weighted average of one-dimensional array
  - Supports missing value handling
  - Used for statistical analysis, data averaging, weighted data processing

## 42. dimsumwgt.f [📝 src](fortran/dimsumwgt.f)
**Function**: One-dimensional weighted sum calculation

**Detailed Description**:
- Contains one subroutine: DIMSUMWGT
- Input parameters: Array x, weights wgt, control parameter iopt
- Output parameters: Weighted sum sumx
- Key features:
  - Calculates weighted sum of one-dimensional array
  - Supports missing value handling
  - Used for statistical analysis, data summation, weighted data aggregation

## 43. distribfit.f [📝 src](fortran/distribfit.f)
**Function**: Weibull distribution parameter fitting

**Detailed Description**:
- Contains two subroutines: WEIBFIT, QNORMALQ
- Input parameters: Data array xin, minimum sample size nmin, confidence level level
- Output parameters: Weibull distribution parameters weib(6)
- Key features:
  - Uses maximum likelihood method to fit two-parameter Weibull distribution
  - Used for reliability analysis, survival analysis, extreme value statistics, failure analysis

## 44. dmapgci.f [📝 src](fortran/dmapgci.f)
**Function**: Great circle interpolation on Earth's surface

**Detailed Description**:
- Contains two subroutines: DMAPGCI, DGCDIST
- Input parameters: Longitude and latitude of two points, number of interpolation points
- Output parameters: Interpolation point coordinates on great circle path
- Key features:
  - Interpolates along great circle path between two points on Earth's surface
  - Calculates great circle distances
  - Used for geographic information systems, navigation route planning, geographic data interpolation

## 45. dpsort.f [📝 src](fortran/dpsort.f)
**Function**: Double precision passive sorting

**Detailed Description**:
- Contains two subroutines: DPSORTDRIVER, DPSORT
- Input parameters: Array dx, sorting control parameter kflag
- Output parameters: Sorting index array iperm
- Key features:
  - Passive sorting algorithm using modified quicksort
  - Generates permutation vector without changing original array
  - Used for data sorting, index generation, data reordering

## 46. dpsort_large.f [📝 src](fortran/dpsort_large.f)
**Function**: Large array passive sorting

**Detailed Description**:
- Contains two subroutines: DPSORTLARGEDRIVER, DPSORTLARGE
- Input parameters: Large array dx, sorting control parameter kflag
- Output parameters: Sorting index array iperm (64-bit integer)
- Key features:
  - Passive sorting algorithm for large arrays
  - Supports 64-bit indexing
  - Used for big data sorting, ultra-large array processing, high-performance computing

## 47. dtrend_dp.f [📝 src](fortran/dtrend_dp.f)
**Function**: Time series trend removal

**Detailed Description**:
- Contains two subroutines: DDTRNDX, DDTRNDMSG
- Input parameters: Time series data, detrending options
- Output parameters: Detrended data, trend coefficients
- Key features:
  - Time series trend removal
  - Supports linear and quadratic trends
  - Handles missing values, used for time series analysis, climate data processing, trend analysis

## 48. dtrend_lsq_msg_dp.f [📝 src](fortran/dtrend_lsq_msg_dp.f)
**Function**: Least squares quadratic trend removal

**Detailed Description**:
- Contains one subroutine: QDTRNDMSG
- Input parameters: Time series Y, missing value flag, control options
- Output parameters: Quadratic detrended series, trend coefficients
- Key features:
  - Uses least squares method for quadratic trend removal
  - Supports missing values
  - Used for time series preprocessing, nonlinear trend analysis, data denoising

## 49. dz_height.f [📝 src](fortran/dz_height.f)
**Function**: Atmospheric layer thickness calculation

**Detailed Description**:
- Contains three subroutines: DZHGTDRV, DZHEIGHT, DZMONO
- Input parameters: Height field z, surface height zsfc, top height ztop
- Output parameters: Layer thickness dz
- Key features:
  - Calculates atmospheric layer thickness
  - Based on given height field and boundary conditions
  - Used for atmospheric models, vertical coordinate processing, layer thickness calculations

## 50. eof2data.f [📝 src](fortran/eof2data.f)
**Function**: EOF analysis data reconstruction

**Detailed Description**:
- Contains one subroutine:
  - DEOF2DATA: Reconstruct data field from EOF modes and time coefficients
- Input parameters:
  - NEVAL,NPTS,NTIM: Number of modes, spatial points, time points
  - EV: EOF eigenvectors (NPTS×NEVAL)
  - EVTS: Principal component time coefficients (NTIM×NEVAL)
  - XMSG: Missing value code
- Output parameters:
  - X: Reconstructed spatiotemporal data field (NTIM×NPTS)
- Key features:
  - Reconstruction formula: $X(t,s) = \sum_k[\text{EVTS}(t,k) \times \text{EV}(s,k)]$
  - Supports partial mode reconstruction (using first NEVAL modes)
  - Automatically handles missing value situations
  - Suitable for data recovery and reconstruction verification after EOF dimensionality reduction

## 51. eof_scripps.f [📝 src](fortran/eof_scripps.f)
**Function**: Empirical Orthogonal Function (EOF) analysis

**Detailed Description**:
- Contains multiple subroutines:
  - DEOF11: Compute EOF analysis for a single dataset using covariance matrix
  - DEOF22: Compute cross-EOF analysis between two datasets
  - DEOF: Main EOF computation program
  - DEOFCOVCOR: Covariance/correlation matrix calculation support
  - DSYMMEIG: Symmetric matrix eigenvalue decomposition
- Input parameters:
  - X: Spatiotemporal data field (time×space)
  - NTIME,NSPACE: Time and space dimensions
  - NEVAL: Number of eigenvectors to compute
  - XMSG: Missing value code
  - IOPT: Calculation option (covariance=0, correlation matrix=1)
- Output parameters:
  - EVAL: Eigenvalues (variance contribution)
  - EVEC: Eigenvectors (EOF modes)
  - PCF: Principal component time coefficients
  - PCVAR: Variance contribution percentage for each mode
- Key features:
  - Based on classical EOF decomposition theory
  - Supports missing value handling
  - Provides both covariance and correlation matrix analysis modes
  - Calculates variance contribution and principal components
  - Suitable for climate mode analysis, dimensionality reduction, pattern recognition

## 52. eqthecalc.f90 [📝 src](fortran/eqthecalc.f90)
**Function**: Equivalent Potential Temperature calculation

**Detailed Description**:
- Contains one main program:
  - DEQTHECALC: WRF model equivalent potential temperature calculation program
- Input parameters:
  - QVP: Water vapor mixing ratio (g/kg), 3D array (miy×mjx×mkzh)
  - TMK: Temperature (K), 3D array (miy×mjx×mkzh)
  - PRS: Total pressure (hPa), 3D array (miy×mjx×mkzh)
  - MIY,MJX,MKZH: Array dimensions
- Output parameters:
  - ETH: Equivalent potential temperature (K), 3D array (miy×mjx×mkzh)
- Key features:
  - Uses WRF physical constants module (EPS, GAMMA, etc.)
  - Calculates lifting condensation level TLCL
  - Equivalent potential temperature formula: $\theta_e = T \times \left(\frac{1000}{p}\right)^{\gamma} \times \exp[\text{thermodynamic term}]$
  - Uses OpenMP parallelization (COLLAPSE(3))
  - Suitable for atmospheric convection analysis, instability diagnostics

## 53. erf_sub.f [📝 src](fortran/erf_sub.f)
**Function**: Error function calculation simple interface

**Detailed Description**:
- Contains two subroutines:
  - DERRF: Wrapper calling standard erf function
  - DERRCF: Wrapper calling erfc1 function
- Input parameters:
  - DERRF: X (scalar argument)
  - DERRCF: IND (integer parameter), X (argument)
- Output parameters:
  - RESULT: Error function calculation result
- Key features:
  - Provides simple NCL interface to standard error functions
  - $\text{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}dt$
  - $\text{erfc}(x) = 1 - \text{erf}(x)$
  - Single-point calculation rather than vectorized
  - Used for probability calculations, statistical analysis

## 54. evssm_lapack.f [📝 src](fortran/evssm_lapack.f)
**Function**: Symmetric matrix eigenvalue/eigenvector calculation

**Detailed Description**:
- Contains one subroutine:
  - EVSSM: Eigenvalue solver using LAPACK DSPEVX
- Input parameters:
  - NVAR: Matrix order
  - LSSM: Symmetric storage matrix size
  - CSSM: Symmetric matrix (packed storage)
  - NEVAL: Number of eigenvalues to compute
  - IOPT: Sorting option (0=descending, 1=ascending)
- Output parameters:
  - WEVAL: Eigenvalue array
  - WEVEC: Eigenvector matrix (NVAR×NEVAL)
  - INFO: LAPACK return status
- Key features:
  - Calls LAPACK DSPEVX routine
  - Computes largest NEVAL eigenvalues and corresponding eigenvectors
  - Supports ascending/descending ordering selection
  - Suitable for principal component analysis, covariance matrix decomposition

## 55. exptapersh.f [📝 src](fortran/exptapersh.f)
**Function**: Spherical harmonic coefficient exponential taper filtering

**Detailed Description**:
- Contains two subroutines:
  - DEXPTAPER: Calculate exponential taper filtering weights
  - DEXPTAPERSH: Apply filtering to spherical harmonic coefficients
- Input parameters:
  - DEXPTAPER: N0 (cutoff parameter), R (exponent), NLAT (number of latitudes)
  - DEXPTAPERSH: AB,BB (spherical harmonic coefficients), N0, R
- Output parameters:
  - DEXPTAPER: S (filtering weight array)
  - DEXPTAPERSH: Filtered AB,BB coefficients
- Key features:
  - Filter weights: $S(n) = \exp\left(-\left(\frac{n(n+1)}{N_0(N_0+1)}\right)^R\right)$
  - Apply weights S(n) to spherical harmonic coefficients AB(m,n) and BB(m,n)
  - Suppress high-order spherical harmonics, preserve low-order components
  - Used for spectral smoothing and denoising of global data

## 56. ezfft_dp.f [📝 src](fortran/ezfft_dp.f)
**Function**: FFTPACK double precision FFT subset

**Detailed Description**:
- Contains multiple subroutines:
  - DEZFFTF: Simplified real forward FFT (outputs Fourier coefficients)
  - DEZFFTB: Simplified real inverse FFT (reconstructs from Fourier coefficients)
  - DEZFFTI: EZ-FFT initialization routine
  - DRFFTF/DRFFTB: Standard real FFT forward/inverse transforms
  - DRFFTI: Real FFT initialization
- Input parameters:
  - N: Sequence length
  - R: Input real sequence
  - WSAVE: Pre-computed work array
- Output parameters:
  - AZERO: DC component
  - A,B: Cosine and sine coefficient arrays
- Key features:
  - Double precision version of FFTPACK library
  - Supports mixed-radix FFT algorithm (2,3,4,5 factor decomposition)
  - EZFFT provides user-friendly coefficient format
  - Efficient real FFT specialized algorithm
  - Suitable for spectral analysis of real time series

## 57. f2foshv.f [📝 src](fortran/f2foshv.f)
**Function**: Vector field transformation from fixed grid to fixed offset grid

**Detailed Description**:
- Contains one main routine:
  - DF2FOSHV: Calls SPHEREPACK's VSHIFTI and VSHIFTE routines
- Input parameters:
  - UOFF,VOFF: u,v components on fixed offset grid (ILON×JLAT)
  - UREG,VREG: u,v components on fixed grid (ILON×JLAT1)
  - ILON,JLAT: Grid dimensions
  - WORK,WSAVE: Work arrays
- Output parameters:
  - IER: Error status code
- Key features:
  - Based on SPHEREPACK library spherical vector field transformation
  - Handles correct transformation of vector components in spherical coordinates
  - Maintains physical consistency of vector fields
  - Used for atmospheric and oceanic model grid transformations

## 58. fill_msg_grid.f [📝 src](fortran/fill_msg_grid.f)
**Function**: Poisson equation method for filling grid missing values

**Detailed Description**:
- Contains two subroutines:
  - POISXY1: Entry routine for missing value filling
  - POISXY2: Poisson equation relaxation iterative solver
- Input parameters:
  - XIO: 2D array with missing values (MX×NY)
  - XMSG: Missing value identifier
  - GTYPE: Boundary condition (0=non-cyclic, 1=x-direction cyclic)
  - NSCAN: Maximum iteration count
  - EPSX,RELC: Convergence criterion and relaxation constant
  - GUESS: Initial guess (0=zero values, 1=zonal average)
- Output parameters:
  - Filled XIO array, MSCAN: Actual iteration count
- Key features:
  - Uses Laplace equation: $\nabla^2 Z = 0$
  - Four-point average relaxation iteration method
  - Supports cyclic and non-cyclic boundary conditions
  - Suitable for spatial data interpolation and missing value imputation

## 59. filtrx.f [📝 src](fortran/filtrx.f)
**Function**: Lanczos and Gaussian digital filter design

**Detailed Description**:
- Contains three subroutines:
  - DFILTRQ: Lanczos filter weight and response function calculation
  - DFILWTQ: Filter weight calculation and normalization
  - FILWGTNORMAL: Gaussian (normal) filter weight generation
- Input parameters:
  - DFILTRQ: NWT(weight count), FCA/FCB(cutoff frequencies), IHP(filter type)
  - FILWGTNORMAL: NWT, SIGMA(standard deviation), ITYPE
- Output parameters:
  - WT: Normalized filter weights
  - FREQ,RESP: Frequency and response arrays
- Key features:
  - Supports low-pass (IHP=0), high-pass (IHP=1), band-pass (IHP=2) filtering
  - Lanczos window reduces Gibbs phenomenon
  - Gaussian filter based on normal distribution
  - Provides complete frequency response analysis
  - Suitable for time series filtering and signal processing

## 60. finfo.f [📝 src](fortran/finfo.f)
**Function**: Fourier series information analysis

**Detailed Description**:
- Contains two subroutines:
  - DFOURINFO: NCL interface for Fourier information analysis
  - DFINFO: Core Fourier series analysis routine
- Input parameters:
  - X: Input time series (NPTS points)
  - NPTS: Number of data points
  - NHRET: Number of harmonics to return
  - SCLPHA: Phase scaling factor
- Output parameters:
  - AMP: Harmonic amplitudes
  - PHA: Harmonic phases (scaled by SCLPHA)
  - PCV: Variance contribution percentage for each harmonic
- Key features:
  - Calls DEZFFTF for Fourier decomposition
  - Amplitude: AMP(n) = $\sqrt{A(n)^2 + B(n)^2}$
  - Phase: PHASE(n) = atan2(B(n),A(n))×adjustment factor
  - Variance contribution: PCV(n) = $(0.5 \times \text{AMP}(n)^2/\text{total variance}) \times 100\%$
  - Suitable for periodic signal analysis and harmonic detection

## 61. fluxEddyTav_dp.f [📝 src](fortran/fluxEddyTav_dp.f)
**Function**: Eddy flux calculation (double precision function)

**Detailed Description**:
- Contains one function:
  - DFLXEDY: Calculate eddy flux between two time series
- Input parameters:
  - X,Y: Two time series arrays (NPTS)
  - NPTS: Number of data points
  - XMSG: Missing value identifier
- Output parameters:
  - DFLXEDY: Returns eddy flux value x'y'
  - IER: Error status code
- Key features:
  - Reynolds decomposition: $x = \bar{X} + x'$, $y = \bar{Y} + y'$
  - Eddy flux: x'y' = x̄ȳ - X̄Ȳ
  - Internal correlation coefficient calculation: $CC = \frac{x'y'}{\sqrt{x'x' \times y'y'}}$
  - Automatic missing value handling
  - Suitable for atmospheric turbulence analysis, heat flux, momentum flux calculations

## 62. fo2fs_dp.f [📝 src](fortran/fo2fs_dp.f)
**Function**: Scalar field transformation from fixed offset grid to fixed grid

**Detailed Description**:
- Contains one main routine:
  - DFO2F: Calls SPHEREPACK's DSSHIFTI and DSSHIFTE routines
- Input parameters:
  - GOFF: Scalar field on fixed offset grid (ILON×JLAT)
  - GREG: Scalar field on fixed grid (ILON×JLAT1)
  - ILON,JLAT: Grid dimensions
  - WORK,WSAVE: Work arrays
- Output parameters:
  - IER: Error status code
- Key features:
  - Based on SPHEREPACK library spherical scalar field transformation
  - Handles correct transformation of scalar fields in spherical coordinates
  - Used for scalar field conversion between different grids
  - Atmospheric and oceanic model data post-processing

## 63. gamfitd3.f [📝 src](fortran/gamfitd3.f)
**Function**: Incomplete gamma distribution parameter estimation

**Detailed Description**:
- Contains one main routine:
  - GAMFITD3: Maximum likelihood method for estimating gamma distribution parameters
- Input parameters:
  - DATARR: Data array (N)
  - N: Number of data points
  - DMSG: Missing value identifier
  - PCRIT: Critical percentage for valid data
  - INVSCALE: Whether to return 1/$\beta$ (0=no, 1=yes)
- Output parameters:
  - ALPHA: $\alpha$ parameter (log-space parameter)
  - SCALE: $\beta$ parameter (scale parameter)
  - SHAPE: $\eta$ parameter (shape parameter)
  - PZERO: Zero value probability
  - IER: Status code
- Key features:
  - Uses Thom (1958) maximum likelihood estimation formula
  - $\alpha = \log(\bar{x}) - (\sum\log(x))/n$
  - $\eta = (1+\sqrt{1+4\alpha/3})/(4\alpha)$
  - Handles zero values and missing values
  - Suitable for precipitation distribution fitting and extreme value analysis

## 64. gamma_interface.f [📝 src](fortran/gamma_interface.f)
**Function**: SLATEC gamma function vectorized interface

**Detailed Description**:
- Contains one main routine:
  - GAMMACOMPLETE: Calls SLATEC DGAMMASLATEC function
- Input parameters:
  - NX: Array length
  - XIN: Input argument array
  - HAS_MSG: Whether missing values exist (0=no, 1=yes)
  - XMSG: Missing value identifier
- Output parameters:
  - XOUT: Gamma function value array
- Key features:
  - Batch calculation of gamma function values
  - Supports missing value handling
  - Calls high-precision SLATEC mathematical library
  - $\Gamma(x) = \int_0^\infty t^{x-1}e^{-t}dt$
  - Suitable for probability distribution calculations and statistical analysis

## 65. gaus.f [📝 src](fortran/gaus.f)
**Function**: SPHEREPACK Gaussian quadrature point calculation

**Detailed Description**:
- Contains four subroutines:
  - GAQDNCL: Main control routine, calculates Gaussian latitudes and weights
  - GAQDNCL1: Sets up tridiagonal matrix and solves eigenvalue problem
  - DRSTNCL: Real symmetric tridiagonal matrix eigenvalue calculation (EISPACK modified version)
  - DINTQLNCL: QL algorithm eigenvalue iterative solver
  - DPYTHANCL: Calculates $\sqrt{a^2+b^2}$ avoiding overflow
- Input parameters:
  - NLAT: Number of Gaussian latitude points
  - WORK: Work array (length $\geq$ 4×NLAT×(NLAT+1)+2)
  - LWORK: Work array length
- Output parameters:
  - THETA: Gaussian latitude points (radians, 0 to $\pi$)
  - WTS: Corresponding Gaussian weights
  - IERROR: Error code
- Key features:
  - Based on Legendre polynomial recurrence relations to construct tridiagonal matrix
  - Eigenvalues are Gaussian points, first eigenvector component squared ×2 gives weights
  - Uses QL algorithm to solve symmetric tridiagonal matrix eigenvalue problem
  - Results sorted in ascending order by radians
  - Used for high-precision numerical integration in spherical harmonic analysis

## 66. gausLobat.f [📝 src](fortran/gausLobat.f)
**Function**: Gauss-Lobatto collocation points and weights calculation

**Detailed Description**:
- Contains multiple subroutines:
  - GAUSLOBAT: NCL interface, returns Gauss-Lobatto latitudes in [-90,90] degrees
  - FINDGLW: Calculates corresponding Gauss-Lobatto weights for given latitudes
  - GAUSSLOBATTO: Core calculation program, solves for Gauss-Lobatto points and weights
  - NEWTONRAPHSON1/2: Hybrid Newton-Raphson and bisection iterative root finding
  - JACOBF: Computes Legendre polynomials and their first and second derivatives
- Input parameters:
  - NPTS: Number of collocation points
  - XGLAT: Input latitude array (for FINDGLW)
- Output parameters:
  - XGL: Gauss-Lobatto collocation points (degrees or [-1,1])
  - WEIGHT: Corresponding weights
- Key features:
  - Based on zeros of Legendre polynomial $P_N(x)$ and its derivative $P'_N(x)$
  - Endpoints included in configuration: $x_1 = -1$, $x_{N+1} = 1$
  - Weight formula: $w_i = \dfrac{2}{N(N+1)\left[P_N(x_i)\right]^2}$
  - High-precision Newton-Raphson iterative solution (tolerance $1\times10^{-15}$)
  - Used for spectral methods and high-precision numerical integration

## 67. gendat.f [📝 src](fortran/gendat.f)
**Function**: Two-dimensional test data generator

**Detailed Description**:
- Contains two subroutines:
  - DGENDAT: Generates 2D test data field with specified characteristics
  - DFRAN: Simple random number generator
- Input parameters:
  - DATA: Output data array (IDIM×N)
  - M,N: Actual array dimensions used
  - MLOW,MHGH: Number of low and high value centers (1-25)
  - DLOW,DHGH: Minimum and maximum data values
  - ISEED: Random number seed (0-100)
- Output parameters:
  - Filled DATA array
- Key features:
  - Uses exponential function to generate data field: $\sum \exp(-\text{distance}^2)$
  - Randomly distributed high and low value centers
  - Precise control of minimum and maximum values
  - 100 preset random number sequences
  - Suitable for graphics program testing and algorithm validation

## 68. hydro_dp.f [📝 src](fortran/hydro_dp.f)
**Function**: Hydrostatic height calculation

**Detailed Description**:
- Contains one subroutine:
  - DHYDRO: Calculates geopotential height using hydrostatic equation
- Input parameters:
  - P: Pressure array (Pa, NLVL levels)
  - TKV: Virtual temperature array (K, NLVL levels)
  - ZSFC: Surface geopotential height (gpm)
  - NLVL: Number of vertical levels
- Output parameters:
  - ZH: Geopotential height at each level (gpm)
  - IER: Error code
- Key features:
  - Hydrostatic equation: dZ = -(RT/g)d(ln p)
  - Integral form: $Z(k) = Z(k-1) + \frac{R}{g} \times \bar{T} \times \ln\left(\frac{p(k-1)}{p(k)}\right)$
  - T̄ uses logarithmic weighted average: T̄ = (T₁ln p₁ + T₂ln p₂)/ln(p₁p₂)
  - Uses WMO standard gravitational acceleration: $g = 9.80665$ m/s$^2$
  - Dry air gas constant: $R = 287.04$ J/(kg·K)
  - Suitable for atmospheric physics height calculations

## 69. hyi2hyo.f [📝 src](fortran/hyi2hyo.f)
**Function**: Hybrid coordinate system vertical interpolation

**Detailed Description**:
- Contains one subroutine:
  - DHYI2HYOB: Interpolates variables between different hybrid coordinate levels
- Input parameters:
  - P0: Reference pressure (Pa)
  - HYAI,HYBI: Input hybrid coordinate coefficients (KLEVI levels)
  - HYAO,HYBO: Output hybrid coordinate coefficients (KLEVO levels)
  - PSFC: Surface pressure (MLON×NLAT)
  - XI: Input variable (MLON×NLAT×KLEVI)
  - INTFLG: Extrapolation flag (0=set missing value, 1=use nearest value)
  - XMSG: Missing value identifier
- Output parameters:
  - XO: Output variable (MLON×NLAT×KLEVO)
  - MSGFLG: Missing value presence flag
- Key features:
  - Hybrid coordinate pressure: $p(k) = \text{hya}(k) \times p_0 + \text{hyb}(k) \times p_{\text{sfc}}$
  - Log-linear interpolation: linear interpolation in ln p space
  - Handles out-of-range cases (extrapolation or set missing values)
  - Suitable for vertical coordinate conversion between different atmospheric models

## 70. index77_dp.f [📝 src](fortran/index77_dp.f)
**Function**: Climate index calculation (such as Southern Oscillation Index)

**Detailed Description**:
- Contains two subroutines:
  - DINDX77: Simplified interface, manages work array allocation
  - DINDX77X: Core calculation program
- Input parameters:
  - X,Y: Monthly data from two stations (NMOS×NYRS)
  - NMOS,NYRS: Number of months and years (typically NMOS=12)
  - XMSG: Missing value identifier
  - IPRNT: Print flag
- Output parameters:
  - ZI: Climate index (NMOS×NYRS)
  - ZNI: Noise index (NMOS×NYRS)
- Key features:
  - Based on Trenberth(1984) method for calculating Southern Oscillation Index
  - Calculates monthly climate means and standard deviations
  - Overall anomaly standard deviation normalization: $\text{ZI} = \frac{X'}{\sigma_X} - \frac{Y'}{\sigma_Y}$
  - Noise index: $\text{ZNI} = \frac{X'}{\sigma_X} + \frac{Y'}{\sigma_Y}$
  - Suitable for ENSO, AO and other climate index calculations

## 71. int2p_dp.f [📝 src](fortran/int2p_dp.f)
**Function**: Variable interpolation between different pressure coordinate levels

**Detailed Description**:
- Contains one subroutine:
  - DINT2P: Pressure coordinate double precision interpolation program
- Input parameters:
  - PPIN,XXIN: Input pressure and variable arrays (NPIN)
  - PPOUT: Target output pressure array (NPOUT)
  - LINLOG: Interpolation type (1=linear, 0=log-linear, negative=allow extrapolation)
  - XMSG: Missing value identifier
- Output parameters:
  - XXOUT: Interpolated variable array
  - IER: Error status code
- Key features:
  - Supports automatic identification and reordering of increasing and decreasing pressure sequences
  - Linear interpolation: $\text{SLOPE} = \frac{X_1-X_2}{P_1-P_2}$
  - Logarithmic interpolation: $\text{SLOPE} = \frac{X_1-X_2}{\ln(P_1)-\ln(P_2)}$
  - Extrapolation option: uses endpoint slopes for boundary extrapolation
  - Automatically handles missing values and exactly matching pressure levels
  - Suitable for atmospheric model vertical coordinate conversion

## 72. julGreg.f [📝 src](fortran/julGreg.f)
**Function**: Mutual conversion between Julian day numbers and Gregorian calendar dates

**Detailed Description**:
- Contains four functions:
  - GREG2JULI: Gregorian calendar to integer Julian day number
  - GREG2JULD: Gregorian calendar to double precision Julian day number (including hours)
  - JULD2GREG: Double precision Julian day number to Gregorian calendar
  - JULI2GREG: Integer Julian day number to Gregorian calendar
- Input parameters:
  - GREG2JULI/JULD: YYYY,MM,DD[,HR] (year, month, day[, hour])
  - JULD2GREG/JULI2GREG: JULD/JULD (Julian day number)
- Output parameters:
  - Gregorian calendar: YYYY,MM,DD,HR
  - Julian day number: INTEGER or REAL*8
- Key features:
  - Based on Fliegel & Van Flandern (1968) algorithm
  - Julian day counting starts from January 1, 4713 BC
  - Julian day starts at noon 12:00 UT
  - Supports hour precision time conversion
  - Suitable for historical and astronomical time calculations

## 73. kens_trend_dp.ncl.f [📝 src](fortran/kens_trend_dp.ncl.f)
**Function**: Kendall's $\tau$ nonparametric trend test and Theil-Sen slope estimation

**Detailed Description**:
- Contains two subroutines:
  - KENSTSTNCL: NCL interface program, handles work array allocation
  - KENSTST: Core Kendall trend test implementation
- Input parameters:
  - XDATA: Time series data (N)
  - N: Number of data points (must be $\geq$10)
  - EPS: Equality judgment threshold (default 1d-5)
- Output parameters:
  - S: Kendall's S statistic
  - Z: Standardized Z value
  - PROB: Significance probability
  - SLOPE: Theil-Sen slope array (for median slope estimation)
- Key features:
  - S statistic: $S = \sum_{j>i} \text{sign}(x_j - x_i)$
  - Z statistic: $Z = \frac{S \pm 1}{\sqrt{\text{Var}(S)}}$, variance adjusted for tied values
  - Significance: $\text{PROB} = \text{erf}\left(\frac{|Z|}{\sqrt{2}}\right)$
  - Theil-Sen slope: median of all pairwise slopes
  - Suitable for trend analysis of non-normal distributions and data with outliers

## 74. kernel_density.f [📝 src](fortran/kernel_density.f)
**Function**: Gaussian kernel density estimation and plug-in bandwidth selection algorithm

**Detailed Description**:
- Contains two subroutines:
  - KERDENI: NCL interface program, handles missing values and data preprocessing
  - PLUGIN: Core plug-in bandwidth selection and density estimation algorithm
- Input parameters:
  - X: Input data array (N)
  - Z: Estimation grid point array (M)
  - XMSG: Missing value identifier
- Output parameters:
  - F: Density estimation values at each grid point (M)
  - H: Optimized bandwidth value
- Key features:
  - Gaussian kernel function: K(u) = exp(-u$^2$/2)/$\sqrt{2\pi}$
  - Plug-in bandwidth selection: iterative optimization for optimal bandwidth
  - Density estimation: $f(x) = \frac{1}{nh}\sum K\left(\frac{x-X_i}{h}\right)$
  - IQR initial bandwidth: $h_0 = 0.92 \times \text{IQR} \times n^{-1/7}$
  - 5 iterations of optimization to obtain final bandwidth
  - Suitable for nonparametric density estimation of continuous distributions

## 75. kmeans_kmns_as136.f [📝 src](fortran/kmeans_kmns_as136.f)
**Function**: K-means clustering algorithm (AS-136 standard implementation)

**Detailed Description**:
- Contains four subroutines:
  - KMNS136: NCL interface program, handles initialization and cluster center setup
  - KMNS: K-means algorithm main program
  - OPTRA: Optimal transfer stage implementation
  - QTRAN: Quick transfer stage implementation
- Input parameters:
  - DAT: Data matrix (M×N), M points with N-dimensional features
  - K: Number of clusters
  - ITER: Maximum number of iterations
  - ISEED: Initialization seed (1=sequential, 2=interval sampling)
- Output parameters:
  - IC1: Cluster assignment for each point (M)
  - CLCNTR: Final cluster centers (N×K)
  - WSS: Within-cluster sum of squares for each cluster (K)
  - IER: Error code (0=success, 1=empty cluster, 2=no convergence, 3=parameter error)
- Key features:
  - Based on Hartigan & Wong (1979) AS-136 algorithm
  - Two-stage iteration: optimal transfer (global optimization) + quick transfer (local optimization)
  - Objective function: minimize within-cluster sum of squares $\sum\sum\|x_i - c_k\|^2$
  - Automatic detection of empty clusters and convergence
  - Suitable for unsupervised clustering analysis of numerical data


## 76. kron_square.f [📝 src](fortran/kron_square.f)
**Function**: Calculation of Kronecker product of two square matrices

**Detailed Description**:
- Contains one subroutine:
  - KRONSQM: Square matrix Kronecker product calculation program
- Input parameters:
  - N: Input matrix dimension (N×N)
  - N2: Output matrix dimension (N$^2$×N$^2$)
  - A,B: Two input square matrices (N×N)
- Output parameters:
  - C: Kronecker product matrix (N$^2$×N$^2$)
- Key features:
  - Kronecker product definition: $C = A \otimes B$
  - Element calculation: $C(\text{row},\text{col}) = A(i,j) \times B(k,l)$
  - Index mapping: $\text{row} = n \times (i-1) + k$, $\text{col} = n \times (j-1) + l$
  - Output matrix dimension is the square of input matrix dimension
  - Suitable for linear algebra, matrix analysis, and system theory

## 77. lclvl.f [📝 src](fortran/lclvl.f)
**Function**: Calculate pressure value corresponding to lifted condensation level (LCL) of air parcel

**Detailed Description**:
- Contains one subroutine:
  - DLCLPRS: Lifted condensation level pressure calculation program
- Input parameters:
  - P: Initial pressure (mb/hPa)
  - TK: Initial temperature (K)
  - TDK: Dew point temperature (K)
- Output parameters:
  - PLCL: Pressure corresponding to lifted condensation level (mb/hPa)
- Key features:
  - Based on Stipanuk (1973) iterative algorithm
  - Saturation vapor pressure: $\text{ES} = \text{ES0} \times \exp\left(\frac{17.67 \times \text{TD}}{\text{TD}+243.5}\right)$
  - Saturation mixing ratio: $\text{WS} = \frac{622 \times \text{ES}}{P - \text{ES}}$  
  - Potential temperature: $\text{PTK} = \text{TK} \times \left(\frac{1000}{P}\right)^{0.286}$
  - 10 iterations to solve intersection of dry adiabat and constant mixing ratio lines
  - Suitable for atmospheric physics, convection analysis, and weather forecasting

## 78. localMinMax.f [📝 src](fortran/localMinMax.f)
**Function**: Detect local minima and maxima in two-dimensional grids

**Detailed Description**:
- Contains two subroutines:
  - DLOCALMN: Local minima detection program
  - DLOCALMX: Local maxima detection program
- Input parameters:
  - X: 2D data grid (NI×NJ)
  - XMSG: Missing value identifier
  - LWRAP: Boundary wrap flag (1=longitude wrap, 0=no wrap)
  - DELTA: Extrema judgment threshold (typically=0.0)
- Output parameters:
  - XI,YJ: i,j coordinates of extrema points (adjusted to NCL indices)
  - LOCMIN/LOCMAX: Values at extrema points
  - NMIN/NMAX: Number of detected extrema points
- Key features:
  - Nine-point template detection: center point compared with surrounding 8 neighbors
  - Minimum condition: all neighboring points > center value + DELTA
  - Maximum condition: all neighboring points < center value - DELTA
  - Supports longitude wrap boundary handling
  - Skips regions containing missing values
  - Suitable for meteorological field extrema analysis and pattern recognition

## 79. lspoly.f [📝 src](fortran/lspoly.f)
**Function**: Weighted least squares polynomial fitting algorithm

**Detailed Description**:
- Contains two subroutines:
  - DLSPOLY: NCL interface program, handles parameter checking and work arrays
  - DLSPLY2: Core least squares fitting algorithm implementation
- Input parameters:
  - X,Y: Data point coordinates (M points)
  - WGT: Weight array (M)
  - N: Number of polynomial coefficients (degree = N-1)
  - M: Number of data points (must have M > N)
- Output parameters:
  - COEF: Polynomial coefficients (N), COEF(1) is constant term
  - IER: Error code (0=success, 1=M<N, 2=N>5 warning, 3=singular matrix)
- Key features:
  - Normal equation method: $\mathbf{A}^T \mathbf{W} \mathbf{A} \times \text{COEF} = \mathbf{A}^T \mathbf{W} \mathbf{Y}$
  - Gaussian elimination solution: complete pivoting strategy
  - Supports 0-4 degree polynomials ($N \leq 5$)
  - Weighted fitting: weight of 0 means exclude that point
  - Numerical stability: high-degree polynomials recommend other methods
  - Suitable for trend fitting, curve regression, and data smoothing

## 80. mixhum_ptd.f [📝 src](fortran/mixhum_ptd.f)
**Function**: Calculate water vapor mixing ratio or specific humidity from pressure, temperature, and dew point temperature

**Detailed Description**:
- Contains one subroutine:
  - DWMRQ: Water vapor mixing ratio/specific humidity calculation program
- Input parameters:
  - P: Pressure array (Pa), converted to mb/hPa
  - TD: Dew point temperature array (K), converted to °C
  - PMSG,TDMSG: Missing value identifiers for pressure and dew point temperature
  - ISWIT: Calculation option (1=mixing ratio kg/kg, 2=specific humidity kg/kg, negative=return g/kg)
- Output parameters:
  - WMR: Water vapor mixing ratio or specific humidity array
  - WMRMSG: Output missing value identifier
- Key features:
  - Calls DWMRSKEWT function to calculate basic mixing ratio
  - Mixing ratio (kg/kg): $\text{WMR} = \text{DWMRSKEWT}(P_{\text{mb}},TD_C) \times 0.001$
  - Specific humidity conversion: $\text{SH} = \frac{\text{WMR}}{1 + \text{WMR}}$
  - Unit conversion: negative ISWIT returns g/kg units
  - Automatically handles missing values
  - Suitable for atmospheric humidity analysis and water vapor transport calculations  
## 81. mixhum_ptrh.f [📝 src](fortran/mixhum_ptrh.f)
**Function**: Calculate water vapor mixing ratio or specific humidity from pressure, temperature, and relative humidity

**Detailed Description**:
- Contains one subroutine:
  - DMIXHUM1: Water vapor mixing ratio/specific humidity calculation program
- Input parameters:
  - P: Pressure (hPa or mb)
  - TK: Temperature (K)
  - RH: Relative humidity (percentage)
  - ISWIT: Calculation option (1=mixing ratio, 2=specific humidity, negative=kg/kg units)
- Output parameters:
  - QW: Water vapor mixing ratio or specific humidity
- Key features:
  - Saturation vapor pressure: $\text{EST} = 6.11 \times \exp\left(\frac{17.269 \times (\text{TK}-273.15)}{\text{TK}-35.86}\right)$
  - Mixing ratio formula: $\text{QW} = \frac{0.622 \times \text{EST}}{P - 0.378 \times \text{EST}} \times \frac{\text{RH}}{100}$
  - Specific humidity conversion: $\text{SH} = \frac{\text{QW}}{1 + \text{QW}}$
  - Unit control: negative ISWIT returns kg/kg, positive returns g/kg
  - Suitable for atmospheric humidity calculations and dew point analysis

## 82. mlegev_memory.f [📝 src](fortran/mlegev_memory.f)
**Function**: Maximum likelihood estimation of generalized extreme value distribution parameters using AS-215 algorithm

**Detailed Description**:
- Contains two subroutines:
  - DMLEGEVI: NCL interface program, handles missing values and parameter initial value setting
  - DMLEGEV: AS-215 core algorithm implementation
- Input parameters:
  - X: Data array (after removing missing values, $\geq$10 points)
  - XMSG: Missing value identifier
  - MONIT: Monitoring flag (0=no monitoring)
- Output parameters:
  - VALS(6): [location parameter $\xi$, scale parameter $\sigma$, shape parameter $\kappa$, standard errors ×3]
  - IERR: Error code (0=success, 1=N<2, 2=no convergence, 3=out of bounds, 4=step size limit)
- Key features:
  - Based on Hosking (1985) AS-215 algorithm
  - Newton-Raphson and steepest ascent hybrid iteration
  - Automatic parameter initial estimation: $\sigma_0=\text{std}\times\sqrt{6/\pi}$, $\xi_0=\text{mean}-\gamma\times\sigma_0$
  - Maximum 30 iterations, 50 likelihood evaluations
  - Suitable for extreme value statistics, risk analysis, and hydrological frequency analysis

## 83. obsp1_mult_time_dp.f [📝 src](fortran/obsp1_mult_time_dp.f)
**Function**: Iterative improvement-based objective analysis interpolation algorithm

**Detailed Description**:
- Contains two subroutines:
  - DOBJANLX: Driver program, handles data preprocessing and sorting
  - DOBJANL: Core objective analysis algorithm implementation
- Input parameters:
  - PLON,PLAT,PVAL: Observation point longitude, latitude and values (NTIM×NPTS)
  - GLAT,GLON: Output grid longitude and latitude (NLAT×MLON)
  - RSCAN: Scan radius sequence (degrees), in decreasing order
  - SMWGT: Blending weights (NSCAN)
  - XMSG,PMSG: Missing value identifiers
- Output parameters:
  - GRID: Interpolated grid (MLON×NLAT×NTIM)
  - IER: Error code
- Key features:
  - Multi-scale analysis: large → small radius progressive refinement
  - Gaussian weights: $\text{wgt} = \exp(-4 \times (\text{distance}/\text{radius})^2)$
  - Distance calculation: great circle distance, exact spherical geometry
  - Iterative blending: $\text{grid} = \text{smwgt} \times \text{new value} + (1-\text{smwgt}) \times \text{old value}$
  - Suitable for meteorological objective analysis and station data gridding

## 84. ocean.f [📝 src](fortran/ocean.f)
**Function**: Area-weighted smoothing and mixed layer depth calculation for ocean data

**Detailed Description**:
- Contains two subroutines:
  - WGT_AREA_SMOOTH: Area-weighted smoothing algorithm
  - MIXED_LAYER_DEPTH: Mixed layer depth calculation
- WGT_AREA_SMOOTH input parameters:
  - FIELD: Input field (NX×NY×NO)
  - AREA: Area weights (NX×NY)
  - ICYCLIC: Longitude cyclic boundary flag
  - FILL_VALUE: Missing value
- MIXED_LAYER_DEPTH input parameters:
  - FIELD: Density field (NX×NY×NZ)
  - KMT: Bottom topography (NX×NY)
  - HT: Sea surface height (NX×NY)
  - DEPTH: Depth array (NZ)
  - OFFSET: Density discrimination threshold
- Output parameters:
  - FIELD_RET: Smoothed field or mixed layer depth
- Key features:
  - Five-point weighted smoothing: center + east-west-north-south neighbors
  - Area weight normalization: sum of weights = 1
  - Mixed layer definition: depth position of surface density + offset
  - Linear interpolation for depth: $\text{depth} = d_1 + (\rho-\rho_1) \times (d_2-d_1)/(\rho_2-\rho_1)$
  - Suitable for ocean circulation analysis and mixed layer studies

## 85. omcalc_ccm.f [📝 src](fortran/omcalc_ccm.f)
**Function**: Diagnostic calculation of vertical pressure velocity (omega) in CCM climate model

**Detailed Description**:
- Contains one subroutine:
  - OMCALCCCM: Omega calculation program modified based on CAM3.0
- Input parameters:
  - U,V: Zonal and meridional wind speeds (ILON×JLAT×KLEV)
  - D: Divergence field (ILON×JLAT×KLEV)
  - DPSL,DPSM: Longitude and latitude components of surface pressure gradient
  - PMID,PDEL: Mid-layer pressure and layer thickness (ILON×JLAT×KLEV)
  - PSFC: Surface pressure (ILON×JLAT)
  - HYBD,HYBM: Hybrid coordinate coefficients (KLEV)
  - NPRLEV: First pure pressure level
- Output parameters:
  - OMEGA: Vertical pressure velocity (ILON×JLAT×KLEV)
- Key features:
  - First calculate $\omega/p$, then scale to $\omega$
  - Hydrostatic matrix equation: $\text{HKK} \times \omega/p + \text{HLK} \times \sum \text{divergence} \times \Delta p$
  - Hybrid coordinate contribution: handling of $\vec{v} \cdot \nabla p_s$ term
  - Layer-by-layer integration: calculation from top to bottom layers
  - Final scaling: $\omega = (\omega/p) \times p_{\text{mid}}$
  - Suitable for atmospheric dynamics diagnosis and vertical motion analysis  
## 86. p2hyo.f [📝 src](fortran/p2hyo.f)
**Function**: Interpolate constant pressure level data to hybrid coordinate layers

**Detailed Description**:
- Contains two subroutines:
  - P2HYO: Main interpolation interface program, handles parameter checking and array validation
  - P2HYB: Core interpolation algorithm implementation
- Input parameters:
  - PI: Input pressure level array (KLEVI), must be monotonic
  - XI: Input variable (MLON×NLAT×KLEVI)
  - PSFC: Surface pressure (MLON×NLAT), units in Pa
  - P0: Reference pressure, reference pressure for HYAO coefficients
  - HYAO,HYBO: Output hybrid coordinate coefficients (KLEVO)
  - KFLAG: Extrapolation option (0-4, controls out-of-bounds handling)
- Output parameters:
  - XO: Interpolated variable (MLON×NLAT×KLEVO)
  - IFLAG: Missing value presence flag
  - IER: Error code (0=success, other=pressure array not monotonic)
- Key features:
  - Hybrid coordinate pressure: $p(k) = \text{hyao}(k) \times p_0 + \text{hybo}(k) \times p_{\text{sfc}}$
  - Log-linear interpolation: linear interpolation in ln(p) space
  - 5 extrapolation modes: no extrapolation, nearest value, one-sided, two-sided extrapolation
  - Boundary handling: DXDP gradient extrapolation formula
  - Suitable for climate model vertical coordinate conversion

## 87. paleo_coasts.f [📝 src](fortran/paleo_coasts.f)
**Function**: Generate EZMAP paleogeographic coastline database from mask data

**Detailed Description**:
- Contains multiple subroutines:
  - PALEOOUTLINE: Main program, coordinates mask conversion and file generation
  - SVBLED: Core boundary line tracing algorithm
  - WTEMLR: EZMAP line data file writing
  - WRNUMB,WRCHAR: Number and character formatted output
- Input parameters:
  - MASK: Land-sea mask array (NLON×NLAT), MSKVAL=land
  - LAT,LON: Latitude and longitude arrays (NLAT,NLON)
  - ZDAT: Temporary array (IM×JM), IM=2×NLON+1, JM=2×NLAT+1
  - FLINES,FNAMES: Output line data and name data file names
  - MSKVAL: Land mask value
- Output parameters:
  - Generated EZMAP database files
- Key features:
  - Boundary tracing: 8-direction vector search algorithm
  - Coordinate transformation: linear mapping from mask grid to geographic coordinates
  - File format: standard EZMAP binary format
  - Region naming: left side land ("Land"), right side ocean ("Water")
  - Suitable for paleoclimate research and geological history reconstruction

## 88-92. PASSB Series - FFT Backward Pass Function Set
**Function**: Pass functions for FFT backward transform (inverse transform) in FFTPACK

**Detailed Description**:
- PASSB: General backward pass function, handles arbitrary prime factors
- Input parameters:
  - NAC,IDO,IP,L1,IDL1: FFT algorithm control parameters
  - CC,C1,C2,CH,CH2: Input/output arrays
  - WA: Precomputed trigonometric function table
- Output parameters:
  - Transformed complex arrays
- Key features:
  - Mixed radix algorithm: handles prime factor decomposition of 2,3,4,5, etc.
  - Double buffering technique: CC↔CH arrays used alternately
  - Trigonometric function optimization: pre-storage avoids repeated calculations
  - Symmetry utilization: IPPH=(IP+1)/2 half-length calculation
  - Suitable for complex to time-domain FFT inverse transform

## 93-97. PASSF Series - FFT Forward Pass Function Set  
**Function**: Pass functions for FFT forward transform in FFTPACK

**Detailed Description**:
- PASSF: General forward pass function, handles arbitrary prime factors
- Input parameters:
  - NAC,IDO,IP,L1,IDL1: FFT algorithm control parameters
  - CC,C1,C2,CH,CH2: Input/output arrays
  - WA: Precomputed trigonometric function table
- Output parameters:
  - Transformed complex arrays
- Key features:
  - Symmetric to PASSB algorithm: sign differences at key operations
  - Forward transform: time domain to frequency domain conversion
  - Efficient decomposition: large-length FFT decomposed into small-length FFT products
  - Bit reversal: natural order to bit-reversed order conversion
  - Suitable for time domain to complex frequency domain FFT transform

## 98. patternCor.f [📝 src](fortran/patternCor.f)
**Function**: Calculate pattern (anomaly) correlation coefficient between two two-dimensional fields

**Detailed Description**:
- Contains two subroutines:
  - PATCOR1: Pattern correlation calculation using one-dimensional weights
  - PATCOR2: Pattern correlation calculation using two-dimensional weights
- Input parameters:
  - X,Y: Two data fields (MLON×NLAT)
  - W: Weight array (NLAT or MLON×NLAT)
  - XMSG,YMSG: Missing value identifiers
- Output parameters:
  - R: Pattern correlation coefficient
  - IER: Error code (0=success, 1=all missing, 2=constant field)
- Key features:
  - Weighted average: $\bar{x} = \frac{\sum(w \times x)}{\sum w}$, $\bar{y} = \frac{\sum(w \times y)}{\sum w}$
  - Anomaly fields: $x' = x-\bar{x}$, $y' = y-\bar{y}$
  - Correlation coefficient: $r = \frac{\sum(w \times x' \times y')}{\sqrt{\sum(w \times x'^2) \times \sum(w \times y'^2)}}$
  - Weight support: latitude weights (cos φ) or area weights
  - Suitable for climate model evaluation and field similarity analysis

## 99. patternCor3.f [📝 src](fortran/patternCor3.f)
**Function**: Calculate pattern correlation coefficient between three-dimensional data fields

**Detailed Description**:
- Three-dimensional version extended from patternCor.f
- Key features:
  - Supports three-dimensional weight distribution
  - Vertical layer correlation analysis
  - Spatiotemporal pattern recognition
  - Multi-layer atmospheric/oceanic data correlation
  - Suitable for 3D climate field comparison and model validation

## 100. phybrid_ccm.f [📝 src](fortran/phybrid_ccm.f)
**Function**: Pressure field calculation for CCM climate model hybrid coordinate system

**Detailed Description**:
- Based on CCM/CAM hybrid vertical coordinate system
- Key features:
  - Hybrid coordinate formula: $p(i,j,k) = \text{hyam}(k) \times p_0 + \text{hybm}(k) \times p_s(i,j)$
  - Layer interface pressure: $p_{\text{interface}} = \text{hyai}(k) \times p_0 + \text{hybi}(k) \times p_s(i,j)$
  - Layer thickness calculation: $\Delta p = p_{\text{interface}}(k+1) - p_{\text{interface}}(k)$
  - Terrain following: near-surface $\sigma$ coordinates, upper-level pressure coordinates
  - Suitable for atmospheric model diagnosis and vertical coordinate processing  

## 101. plotfmt_rddata.f [📝 src](fortran/plotfmt_rddata.f)
**Function**: Plot format data reading interface

**Detailed Description**:
- Contains one subroutine:
  - PLOTFMT_RDDATA: Read two-dimensional data array from fixed file unit
- Input parameters:
  - NX,NY: Data array dimensions
  - FUNIT: File unit number fixed at 10
- Output parameters:
  - SLAB: Read two-dimensional data array (NX×NY)
  - ISTATUS: Status code (0=success, 1=read error)
- Key features:
  - Uses Fortran unformatted binary reading
  - Simple error handling mechanism
  - Fixed file unit number design
  - Suitable for data input interface of plotting programs

## 102. potmp_dpth2pres.f [📝 src](fortran/potmp_dpth2pres.f)
**Function**: Oceanographic potential temperature calculation and depth-to-pressure conversion

**Detailed Description**:
- Contains two subroutines:
  - DPOTMP: Seawater potential temperature calculation
  - DPTH2PRES: Depth to pressure coordinate conversion
- DPOTMP input parameters:
  - PRESS: Pressure (decibar)
  - TEMP: Temperature (Celsius)
  - S: Salinity (PSS-78 standard)
  - RP: Reference pressure (decibar)
- DPTH2PRES input parameters:
  - DEPTH: Depth (meters)
  - HAS_MSG_DEPTH: Missing value flag
  - DEPTH_MSG: Depth missing value code
- Output parameters:
  - POTEMP: Potential temperature (Celsius)
  - PRESSURE: Pressure (bar)
- Key features:
  - Based on Fofonoff (1976) potential temperature calculation formula
  - Uses Levitus (1994) global average temperature-salinity data
  - 4th-order Runge-Kutta integration method
  - Hydrostatic equilibrium principle: $dp = -\rho g dz$
  - Pressure formula: $p = 0.059808 \times (e^{-0.025z}-1) + 0.100766z + 2.28405 \times 10^{-7} \times z^2$
  - Suitable for oceanographic physical oceanography analysis

## 103. prcwat.f [📝 src](fortran/prcwat.f)
**Function**: Precipitable water (column total water vapor) calculation

**Detailed Description**:
- Contains one subroutine:
  - DPRCWATDP: Calculate precipitable water from specific humidity and layer thickness
- Input parameters:
  - Q: Specific humidity array (kg/kg, KLVL layers)
  - DP: Layer thickness array (Pa, KLVL layers)
  - KLVL: Number of vertical layers
  - QMSG,DPMSG: Missing value codes for specific humidity and layer thickness
- Output parameters:
  - PRCWAT: Precipitable water (kg/m$^2$ or mm)
- Key features:
  - Integration formula: $\text{PRCWAT} = \frac{1}{g}\sum[q(k) \times |dp(k)|]$
  - Uses standard gravitational acceleration: $g = 9.81$ m/s$^2$
  - Automatically handles missing value layers
  - Unit conversion: kg/m$^2$ = mm water column height
  - Suitable for atmospheric water vapor analysis and precipitation potential assessment

## 104. pres_hybrid_jra55_dp.f [📝 src](fortran/pres_hybrid_jra55_dp.f)
**Function**: JRA55 reanalysis hybrid coordinate pressure calculation

**Detailed Description**:
- Contains one subroutine:
  - DPHYBRIDJRA55: JRA55-specific hybrid coordinate pressure calculation
- Input parameters:
  - HYA,HYB: JRA55 interface hybrid coordinate coefficients (KLEVI layers)
  - PSFC: Surface pressure (MLON×NLAT, Pa)
  - MLON,NLAT,KLEV,KLEVI: Grid and layer number dimensions
  - PI: Interface pressure working array
  - PMSG: Missing value code
- Output parameters:
  - PHY: Mid-layer pressure (MLON×NLAT×KLEV, Pa)
  - PI: Interface pressure (MLON×NLAT×KLEVI, Pa)
- Key features:
  - Interface pressure: $p_i(k) = \text{hya}(k) + \text{hyb}(k) \times p_{\text{sfc}}$
  - Mid-layer pressure log-weighted formula: $p = \exp\left[\frac{1}{\Delta p} \times (p_1 \ln p_1 - p_2 \ln p_2) - 1\right]$
  - $\Delta p = p_{\text{interface}}(k) - p_{\text{interface}}(k+1)$
  - Bottom layer forced to 10Pa
  - Based on JRA55 reanalysis system documentation
  - Suitable for vertical coordinate processing of JRA55 data

## 105. preshybrid.f [📝 src](fortran/preshybrid.f)
**Function**: Hybrid coordinate pressure calculation (general version)

**Detailed Description**:
- Contains two subroutines:
  - PRESHYBRID: Hybrid coordinate mid-layer pressure calculation
  - DPRESHYBRID: Hybrid coordinate layer thickness calculation
- PRESHYBRID input parameters:
  - $P_0$: Reference pressure (Pa)
  - $P_s$: Surface pressure (Pa)
  - HYA, HYB: Hybrid coordinate coefficients (KLVL layers)
- DPRESHYBRID input parameters:
  - $P_0$, $P_s$: Reference pressure and surface pressure (Pa)
  - HYAI, HYBI: Interface hybrid coordinate coefficients (KLVL layers)
- Output parameters:
  - PHY: Mid-layer pressure (KLVL, Pa)
  - DPHY: Layer thickness (KLVL-1, Pa)
- Key features:
  - Standard hybrid coordinate formula: $p(k) = \text{hya}(k) \times P_0 + \text{hyb}(k) \times P_s$
  - Layer thickness calculation: $\Delta p(k) = \left| P_0 \times [\text{hyai}(k+1) - \text{hyai}(k)] + [\text{hybi}(k+1) - \text{hybi}(k)] \times P_s \right|$
  - Supports CCM2/CCM3 hybrid coordinate systems
  - Ordered from model top to surface
  - Suitable for atmospheric model vertical coordinate conversion

## 106. prneof_dp.f [📝 src](fortran/prneof_dp.f)
**Function**: Principal component EOF analysis (double precision LAPACK implementation)

**Detailed Description**:
- Contains two subroutines:
  - DDRVEOF: Main EOF analysis driver program
  - DNCLDRV: NCL interface driver program with missing value threshold control
- Input parameters:
  - X: Spatiotemporal data matrix (NOBS×MSTA)
  - NOBS,MSTA: Number of observations and stations
  - XMSG: Missing value code
  - NEVAL: Number of eigenvalues/eigenvectors to calculate
  - JOPT: Matrix type (0=covariance matrix, 1=correlation matrix)
  - PCRIT: Valid data percentage threshold (DNCLDRV)
- Output parameters:
  - EVAL: Eigenvalues (descending order)
  - EVEC: Eigenvectors (spatial modes)
  - PCVAR: Variance contribution percentage for each mode
  - TRACE: Trace of covariance/correlation matrix
- Key features:
  - Based on LAPACK DSPEVX high-precision eigenvalue solver
  - Supports both covariance and correlation matrix analysis
  - Automatically handles missing values and data quality control
  - Eigenvalues sorted in descending order by variance contribution
  - Can calculate principal component time coefficients (requires additional computation)
  - Suitable for climate mode analysis, dimensionality reduction, and pattern recognition

## 107. prneofTranspose.f [📝 src](fortran/prneofTranspose.f)
**Function**: EOF analysis matrix transpose optimization version

**Detailed Description**:
- Contains two subroutines:
  - TDRVPRC: Main interface for transpose EOF analysis
  - XRVEOFT: Core algorithm for transpose matrix EOF
- Input parameters:
  - X: Original data matrix (NROW×NCOL)
  - NROBS,NCSTA: Number of observations and stations
  - XMSG: Missing value code
  - NEVAL: Number of eigenvalues
  - PCRIT: Data completeness threshold percentage
  - JOPT: Standardization option (0=covariance, 1=correlation)
- Output parameters:
  - EVAL: Eigenvalues
  - EVEC: Reconstructed spatial eigenvectors
  - PCVAR: Variance contribution percentage
  - TRACE: Trace of covariance matrix
- Key features:
  - Based on data matrix transpose to improve computational efficiency
  - Significantly accelerates computation when NOBS<<NCSTA
  - Utilizes spatiotemporal duality: EOF_space = f(PC_time, Data_transpose)
  - Automatic data standardization and outlier handling
  - Quality control: removes stations with valid data < PCRIT
  - Automatic work array allocation, no need for preset memory
  - Suitable for large-scale grid data EOF analysis

## 108. psigma.f [📝 src](fortran/psigma.f)
**Function**: $\sigma$ coordinate system pressure calculation

**Detailed Description**:
- Contains one subroutine:
  - DPSIGMA: $\sigma$ coordinate to pressure coordinate conversion
- Input parameters:
  - SIGMA: $\sigma$ coordinate coefficients (KLEV layers)
  - PSFC: Surface pressure (MLON×NLAT, Pa)
  - MLON,NLAT,KLEV: Longitude, latitude, and layer dimensions
- Output parameters:
  - PSIG: $\sigma$ layer pressure (MLON×NLAT×KLEV, Pa)
- Key features:
  - Simple $\sigma$ coordinate formula: $p(i,j,k) = \sigma(k) \times p_{\text{sfc}}(i,j)$
  - $\sigma$ coordinate definition: $\sigma = p/p_{\text{sfc}}$, range [0,1]
  - $\sigma=0$ corresponds to model top, $\sigma=1$ corresponds to surface
  - Terrain-following coordinate system
  - Maintains vertical resolution variation with terrain
  - Suitable for atmospheric models in complex terrain areas
  - Simpler and more direct than hybrid coordinate systems

## 109. random_dp.f [📝 src](fortran/random_dp.f)
**Function**: Double precision random number generation library (complete statistical distributions)

**Detailed Description**:
- Contains multiple random number generation functions:
  - DGENCHI: Chi-square distribution random number generation
  - DGENGAM: Gamma distribution random number generation  
  - DGENNOR: Normal distribution random number generation
  - DGENUNF: Uniform distribution random number generation
  - IGNLGI: Large integer uniform random number generation
  - IGNUIN: Specified range integer random numbers
- Basic random number engine:
  - DRANFD: Basic [0,1) uniform random number generator
  - Based on L'Ecuyer & Cote (1991) algorithm
  - Supports 32 independent random number streams
  - Period length approximately $2.3 \times 10^{18}$
- Management functions:
  - SETALL/SETSD: Set random number seed
  - GETSD: Get current seed state
  - SETCGN/DGETCGN: Switch random number generator
- Key features:
  - High-quality L'Ecuyer linear congruential generator
  - Complete probability distribution function library
  - Supports multi-stream parallel random number generation
  - Statistically validated algorithm implementation
  - Suitable for Monte Carlo simulation and statistical sampling

## 110. rdsstoi.f [📝 src](fortran/rdsstoi.f)
**Function**: NOAA optimal interpolation sea surface temperature data reading

**Detailed Description**:
- Contains two subroutines:
  - RDSSTOI: Read SST data in SSTOI format
  - WRSSTOI: Write Gaussian grid SST data
- RDSSTOI input parameters:
  - NYRSTRT,NYRLAST: Data year range
  - MLON,NLAT: Grid dimensions (fixed 360×180)
- RDSSTOI output parameters:
  - SST: Sea surface temperature array (360×180, °C)
  - INFO: Date and status information (9 integers)
- Data format characteristics:
  - 1°×1° global grid (360 longitude × 180 latitude)
  - Geographic range: 179.5°W-179.5°E, 89.5°S-89.5°N
  - Supports weekly, monthly and climatological data
  - Original data in integer format of °C×100
  - Automatic conversion to floating point temperature values
- Key features:
  - Sequential reading of formatted file 'SSTOI'
  - Each grid contains 8 header information items
  - Header format: (start date, end date, days, index)
  - Data quality identifier: $\leq -1.78$°C indicates sea ice
  - Requires use with land-sea mask
  - Suitable for oceanographic analysis and climate research

## 111. regcoef_dp.f [📝 src](fortran/regcoef_dp.f)
**Function**: Linear regression coefficient calculation and statistical testing

**Detailed Description**:
- Contains one subroutine:
  - DREGCOEF: Complete linear regression analysis
- Input parameters:
  - X,Y: Input data vectors (NPTS points)
  - NPTS: Total number of data points
  - XMSG,YMSG: Missing value codes for X and Y
- Output parameters:
  - RCOEF: Regression coefficient (slope)
  - TVAL: t-statistic (for testing regression coefficient significance)
  - NPTXY: Actual number of data points used
  - XAVE,YAVE: Means of X and Y
  - RSTD: Standard error of regression coefficient
  - YINT: y-intercept
  - IER: Error code
- Key features:
  - Based on Brownlee (1965) statistical theory
  - Regression equation: $Y = \text{RCOEF} \times X + \text{YINT}$
  - t-statistic: $t = \frac{\text{RCOEF}-0}{\text{RSTD}}$, used for significance testing
  - Degrees of freedom: $df = n-2$
  - Automatically handles missing values and boundary conditions
  - Contains complete regression diagnostic statistics
  - Suitable for trend analysis and correlation studies

## 112. relhum_dp.f [📝 src](fortran/relhum_dp.f)
**Function**: Relative humidity calculation (double precision, based on CCM algorithm)

**Detailed Description**:
- Contains one function:
  - DRELHUM: Calculate relative humidity from temperature, mixing ratio, and pressure
- Input parameters:
  - T: Temperature (K)
  - W: Mixing ratio (kg/kg)
  - P: Pressure (Pa)
- Output parameters:
  - DRELHUM: Relative humidity (%)
- Key features:
  - Built-in saturation vapor pressure lookup table (-100°C to +102°C)
  - Temperature range limit: 173.16K to 375.16K
  - Relative humidity formula: $\text{RH} = \frac{W \times (P - 0.378 \times \text{Es})}{0.622 \times \text{Es}} \times 100\%$
  - Linear interpolation for saturation vapor pressure calculation
  - Lower bound constraint: $\text{RH} \geq 0.0001\%$
  - Based on CCM climate model physical processes
  - Suitable for atmospheric humidity analysis and model diagnostics

## 113. relhum_ice.f [📝 src](fortran/relhum_ice.f)
**Function**: Ice surface relative humidity calculation (improved Magnus formula)

**Detailed Description**:
- Contains one subroutine:
  - DRELHUMI: Relative humidity calculation relative to ice surface
- Input parameters:
  - P: Pressure (Pa, internally converted to hPa)
  - TK: Temperature (K)
  - QW: Mixing ratio (kg/kg)
- Output parameters:
  - RH: Relative humidity (%)
- Key features:
  - Based on Alduchov & Eskridge improved Magnus formula
  - Applicable temperature range: +50°C to -80°C
  - Maximum relative error: 0.337%-0.823%
  - Ice surface saturation vapor pressure: $\text{Es} = 6.1128 \times \exp\left(\frac{22.571 \times (T-273.15)}{T-0.44}\right)$
  - Saturation mixing ratio: $\text{Qs} = \frac{0.622 \times \text{Es}}{P - 0.378 \times \text{Es}}$
  - Relative humidity: $\text{RH} = \frac{100 \times \text{QW}}{\text{Qs}}$
  - Suitable for low temperature conditions, upper atmosphere, and polar meteorology

## 114. relhum_water.f [📝 src](fortran/relhum_water.f)
**Function**: Water surface relative humidity calculation (reversible calculation with mixhum_ptrh)

**Detailed Description**:
- Contains one subroutine:
  - DRELHUMW: Relative humidity calculation relative to liquid water
- Input parameters:
  - P: Pressure (Pa, internally converted to hPa)
  - TK: Temperature (K)  
  - QW: Mixing ratio (kg/kg)
- Output parameters:
  - RH: Relative humidity (%)
- Key features:
  - Uses same constants as mixhum_ptrh, ensuring reversible calculation
  - Water surface saturation vapor pressure: $\text{Es} = 6.11 \times \exp\left(\frac{17.269 \times (T-273.15)}{T-35.86}\right)$
  - Saturation mixing ratio: $\text{Qs} = \frac{0.622 \times \text{Es}}{P - 0.378 \times \text{Es}}$
  - Relative humidity: $\text{RH} = \frac{100 \times \text{QW}}{\text{Qs}}$
  - Molecular weight ratio: $M_w/M_d = 0.622$
  - Suitable for normal atmospheric conditions, marine atmosphere, and humid environments

## 115. remap.f [📝 src](fortran/remap.f)
**Function**: Grid remapping with precomputed weights

**Detailed Description**:
- Contains one subroutine:
  - DPOPREMAP: Grid data remapping using precomputed weights
- Input parameters:
  - SRC_ARRAY: Source grid data (NSRC)
  - MAP_WTS: Remapping weight matrix (NW×NLINK)
  - DST_ADD,SRC_ADD: Destination and source grid address indices (NLINK)
  - NLINK: Number of interpolation links
  - NW: Number of weights (typically 1)
  - NDST,NSRC: Number of destination and source grid points
  - XMSG: Missing value code
- Output parameters:
  - DST_ARRAY: Remapped destination grid data (NDST)
- Key features:
  - Based on precomputed sparse weight matrix
  - Weighted interpolation: $\text{DST}(i) = \sum[\text{SRC}(j) \times \text{WGT}(i,j)]$
  - Supports mapping relationships between arbitrary grids
  - Automatically handles missing value propagation
  - Efficient sparse matrix storage format
  - Suitable for conservative interpolation, bilinear interpolation, distance weighting, etc.
  - Commonly used for climate model data post-processing

## 116. rhombtri_dp.f [📝 src](fortran/rhombtri_dp.f)
**Function**: Spherical harmonic coefficient rhombic and triangular truncation (double precision)

**Detailed Description**:
- Contains three subroutines:
  - DRHOMBTRUNC: Rhombic truncation of spherical harmonic coefficients
  - DTRITRUNC: Triangular truncation interface program
  - DTRUNC: Triangular truncation core algorithm
- DRHOMBTRUNC input parameters:
  - AR,AI: Real and imaginary parts of spherical harmonic coefficients (M×N)
  - M,N: Zonal and total wavenumber dimensions
  - R: Rhombic truncation parameter (retained total wavenumber)
- DTRUNC input parameters:
  - A,B: Spherical harmonic coefficient arrays (ID×NM)
  - MS: Truncation parameter
  - ID,NM: Array dimensions
- Output parameters:
  - Truncated spherical harmonic coefficients (coefficients beyond truncation range set to zero)
- Key features:
  - Rhombic truncation: forms rhombic retention region in wavenumber space
  - Triangular truncation: prevents aliasing from product terms
  - Based on SPHEREPACK 2.0 algorithm
  - Supports different cases of $M \geq N$ and $M < N$
  - Used for dealiasing in spherical spectral methods
  - Suitable for spectral transforms in global numerical models

## 117. rmStndx_dp.f [📝 src](fortran/rmStndx_dp.f)
**Function**: Time series mean, median removal and standardization

**Detailed Description**:
- Contains three subroutines:
  - DRMVMEAN: Remove series mean
  - DRMVMED: Remove series median
  - DXSTND: Standardization (zero mean unit variance)
- DRMVMEAN input parameters:
  - X: Input time series (NPTS, modified in place)
  - NPTS: Number of data points
  - XMSG: Missing value code
- DXSTND input parameters:
  - X: Input time series (NPTS, modified in place)
  - IOPT: Standard deviation type (0=sample standard deviation, 1=population standard deviation)
- Output parameters:
  - X: Processed series
  - IER: Error code
- Key features:
  - In-place operation, saves memory
  - Automatically skips missing values
  - Standardization: $X' = (X-\mu)/\sigma$
  - Supports both sample and population standard deviation calculations
  - Suitable for time series preprocessing and anomaly standardization

## 118. rmsd.f [📝 src](fortran/rmsd.f)
**Function**: Root Mean Square Difference (RMSD) calculation

**Detailed Description**:
- Contains one subroutine:
  - DRMSD: Calculate root mean square difference between two arrays
- Input parameters:
  - X,Y: Two input vectors (NPTS)
  - NPTS: Total number of data points
  - XMSG,YMSG: Missing value codes for X and Y
- Output parameters:
  - XYRMSD: Root mean square difference
  - NPTUSE: Actual number of data points used
  - IER: Error code
- Key features:
  - RMSD formula: $\sqrt{\frac{\sum(X-Y)^2}{N}}$
  - Automatically handles missing values
  - Only includes pairs where both arrays have valid values
  - Returns actual number of points used for validation
  - Suitable for model validation, forecast evaluation, and data comparison

## 119. round.f [📝 src](fortran/round.f)
**Function**: Numerical rounding to nearest integer

**Detailed Description**:
- Contains one subroutine:
  - RNDNCL: Round array elements to nearest integer
- Input parameters:
  - XIN: Input floating point array (NPTS)
  - NPTS: Array length
  - ISMSG: Whether missing values exist (0=none, 1=exist)
  - XMSG: Missing value code
  - IOPT: Output option (3=adjust missing values for integer output)
- Output parameters:
  - XOUT: Rounded array (NPTS)
- Key features:
  - Uses ANINT intrinsic function for rounding
  - Automatically handles missing value propagation
  - Automatically adjusts non-integer missing values when IOPT=3
  - Preserves original array (not in-place operation)
  - Suitable for data discretization and integer conversion

## 120. runave_dp.f [📝 src](fortran/runave_dp.f)
**Function**: Running average filtering (double precision)

**Detailed Description**:
- Contains two subroutines:
  - DRUNAVE: Running average main interface
  - DRUNAVX77: Running average core algorithm
- Input parameters:
  - X: Input time series (NPTS, modified in place)
  - NPTS: Number of data points
  - NAVE: Moving window size
  - KOPT: Boundary handling option
  - XMSG: Missing value code
  - WORK: Work array (LWORK)
- Boundary handling options:
  - KOPT<0: Cyclic boundary conditions
  - KOPT=0: Set boundary points to missing values
  - KOPT>0: Reflective (symmetric) boundary conditions
- Output parameters:
  - X: Smoothed time series
  - IER: Error code
- Key features:
  - Supports both odd and even window sizes
  - Three boundary handling methods
  - Automatically handles missing values in the series
  - Any missing value within window causes output point to be missing
  - Suitable for time series smoothing and noise suppression

## 121. runknt.f [📝 src](fortran/runknt.f)
**Function**: Consecutive sequence count statistics

**Detailed Description**:
- Contains two subroutines:
  - NUMRUN2: Two-dimensional array consecutive sequence counting
  - NUMRUN1: One-dimensional array consecutive sequence counting  
- Input parameters:
  - XIN: Input integer array (0/1 values)
  - NX,NTIM: Array dimensions
  - OPTKNT: Counting option (0=overlapping count, 1=unique sequence count)
- Output parameters:
  - XKNT: Count of consecutive sequences of each length (NX×NTIM or NX)
- Key features:
  - OPTKNT=0: Count all consecutive 1-sequences of length n (overlapping)
  - OPTKNT=1: Count independent consecutive 1-sequences (must be bounded by 0s)
  - Automatically pads boundaries with 0s to ensure complete sequence counting
  - Suitable for drought/wet period statistics and consecutive event analysis
  - Climate extreme event duration statistics

## 122. rvdv.f [📝 src](fortran/rvdv.f)
**Function**: Relative vorticity and divergence finite difference calculation

**Detailed Description**:
- Contains four subroutines:
  - DVRFIDF: Relative vorticity finite difference calculation
  - DDVFIDF: Divergence finite difference calculation
  - DLNEXTRP: Corner point linear extrapolation
  - DUSEAVE: Pole point average value processing
- Input parameters:
  - U,V: Wind field components (MLON×NLAT, m/s)
  - GLAT,GLON: Latitude and longitude arrays
  - IOPT: Boundary condition option
  - XMSG: Missing value code
- Output parameters:
  - RV: Relative vorticity (s$^{-1}$)
  - DV: Divergence (s$^{-1}$)
- Key features:
  - Relative vorticity: $\zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + \frac{u}{a}\tan(\phi)$
  - Divergence: $\nabla \cdot \vec{V} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} - \frac{v}{a}\tan(\phi)$  
  - Supports equally-spaced longitude and unequally-spaced latitude
  - Three boundary handling options: cyclic, non-cyclic, mixed
  - Poles use average values, corners use linear extrapolation
  - Suitable for atmospheric dynamics analysis

## 123. s2h.f [📝 src](fortran/s2h.f)
**Function**: $\sigma$ coordinate to hybrid coordinate interpolation

**Detailed Description**:
- Contains two subroutines:
  - DH2SDRV: $\sigma$ to hybrid coordinate conversion driver program
  - DS2HBD: $\sigma$ coordinate interpolation core algorithm
- Input parameters:
  - DATI: $\sigma$ coordinate data (NLVI)
  - HYA,HYB: Hybrid coordinate coefficients (NLVO)
  - P0: Reference pressure (Pa)
  - PSFC: Surface pressure (Pa)
  - SIGI: $\sigma$ coordinate layers (NLVI)
  - INTYP: Interpolation type (1=linear, 2=logarithmic, 3=double logarithmic)
- Output parameters:
  - DATO: Hybrid coordinate data (NLVO)
  - SIGO: Converted $\sigma$ values (NLVO)
- Key features:
  - Hybrid coordinate to $\sigma$: $\sigma(k) = \frac{\text{hya}(k) \times P_0}{P_{\text{sfc}}} + \text{hyb}(k)$
  - Three interpolation methods: linear, logarithmic, double logarithmic
  - Automatically handles model top, bottom, and internal layers
  - Suitable for atmospheric model vertical coordinate transformation

## 124. setrig.f [📝 src](fortran/setrig.f)
**Function**: Lin-Rood model latitude and weight calculation

**Detailed Description**:
- Contains three subroutines:
  - LINROOD: Lin-Rood latitude and weight main interface
  - LINROODWT: Lin-Rood weight calculation
  - SETRIG: Trigonometric function grid setup
- Input parameters:
  - NLAT: Number of latitude points
- Output parameters:
  - LAT: Latitude array (degrees)
  - WEIGHT: Gaussian weight array
  - COSP,COSE,SINP,SINE: Trigonometric function arrays
- Key features:
  - Uniform latitude grid: $\Delta\text{lat} = \frac{180°}{\text{nlat}-1}$
  - Weight calculation: $w(j) = \sin(\phi_{j+1}) - \sin(\phi_j)$
  - Special pole handling: $w(1) = 1 + \sin(\phi_2)$, $w(\text{nlat}) = 1 - \sin(\phi_{\text{nlat}})$
  - Based on NASA Lin-Rood dynamics framework
  - Suitable for finite volume atmospheric models

## 125. sg_tools.f [📝 src](fortran/sg_tools.f)
**Function**: Spherical geometry calculation toolkit

**Detailed Description**:
- Contains multiple spherical geometry functions:
  - GCCONVEX: Determine convexity of spherical polygon
  - GCINOUT: Whether point is inside spherical polygon
  - GCDIST: Shortest distance between two points on sphere
  - GCONARC: Determine if point is on great circle arc
  - GCCONVERT: Spherical angle unit conversion
  - GCINTERP: Great circle interpolation
  - GCQAREA,GCTAREA: Spherical quadrilateral and triangle areas
  - GCPNT2GC: Distance from point to great circle
  - GCDANGLE,GCAANGLE: Angle between great circles
- Input parameters:
  - LAT,LON: Latitude and longitude (degrees)
  - NPTS: Number of points
  - TYPE: Conversion unit type
- Output parameters:
  - Distance, angle, area and other geometric quantities
- Key features:
  - Complete spherical geometry calculation library
  - Supports degrees, radians, meters, kilometers, feet, miles conversion
  - Based on exact spherical trigonometry formulas
  - Automatically handles numerical precision issues
  - Suitable for GIS, geodesy, and spherical grid generation

## 126. shsgc_R42.f [📝 src](fortran/shsgc_R42.f)
**Function**: R42 rhombic truncation spherical harmonic synthesis (hard-coded version)

**Detailed Description**:
- Contains multiple subroutines:
  - DSHSGCR42: R42 spherical harmonic synthesis main interface
  - DTRME: Spherical harmonic coefficient to grid field transformation
  - DALPMN1: Associated Legendre function calculation
  - DCOOL,DFIXRL,DGAUSSV etc.: FFT and Gaussian grid support functions
- Input parameters:
  - A,B: Real and imaginary parts of spherical harmonic coefficients (43×43)
  - WORK: Work array (LWORK)
- Output parameters:
  - GVAR: Grid variable (128×108)
- Key features:
  - Hard-coded R42 rhombic truncation parameters (IRES=42)
  - Fixed 128 longitude × 108 Gaussian latitude grid
  - Based on algorithm provided by Bao Qing from Chinese Academy of Sciences
  - Contains complete FFT transform and Gaussian integration support
  - Uses precomputed Legendre functions for improved efficiency
  - Suitable for T42 global spectral model grid transformation

## 127. simpson_dp.f [📝 src](fortran/simpson_dp.f)
**Function**: Simpson integration rule (double precision, supports missing values)

**Detailed Description**:
- Contains three subroutines:
  - DSIMPEQ: Equal spacing Simpson integration main interface
  - DSIMPNE: Unequal spacing Simpson integration (supports missing values)
  - DSIMPNE2: Unequal spacing Simpson integration core algorithm
- Input parameters:
  - X,Y: Independent and dependent variable arrays (NUM points)
  - DX: Step size for equal spacing case
  - YMSG: Missing value code for Y
- Output parameters:
  - ANS: Integration result
- Key features:
  - Based on analytical integration of 3-point Lagrange interpolation polynomial
  - Equal spacing: error $\propto h^4 f^{(4)}$, where h is step size
  - Unequal spacing: error $\propto (\text{maximum adjacent spacing})^4 \times f^{(4)}$
  - Automatically skips missing value pairs
  - Front-end driver program simplifies interface
  - Based on modified version of NCAR NSSL algorithm

## 128. skewT.f [📝 src](fortran/skewT.f)
**Function**: Skew-T/log-P diagram thermodynamic calculations and plotting tools

**Detailed Description**:
- Contains multiple thermodynamic functions:
  - DSKEWTY: Pressure to Skew-T coordinate Y-axis conversion
  - DSKEWTX: Temperature to Skew-T coordinate X-axis conversion
  - DTMRSKEWT: Mixing ratio line temperature calculation
  - DTDASKEWT: Dry adiabatic line temperature calculation
  - DSATLFTSKEWT: Moist adiabatic line temperature calculation
  - DPTLCLSKEWT: Lifting condensation level calculation
  - DSHOWALSKEWT: Showalter stability index
  - DPWSKEWT: Precipitable water calculation
  - DCAPETHERMO: CAPE convective available potential energy
- Coordinate transformation formulas:
  - $Y = 132.182 - 44.061 \times \log_{10}(P)$
  - $X = 0.54 \times T + 0.90692 \times Y$
- Key features:
  - Complete atmospheric thermodynamic diagnostic tool
  - Based on Stipanuk(1973) and standard meteorological formulas
  - Supports dry and moist adiabatic process calculations
  - Suitable for atmospheric stability analysis and convection forecasting

## 129. slatec_DPOLFT.f [📝 src](fortran/slatec_DPOLFT.f)
**Function**: SLATEC polynomial least squares fitting (double precision)

**Detailed Description**:
- Contains three main subroutines:
  - POLFT/POLFTMSG: NCL interface wrapper programs
  - DPOLFT: SLATEC core polynomial fitting algorithm
  - DP1VLU: Polynomial and derivative evaluation
  - DPCOEF: Orthogonal polynomial coefficient conversion to power series
- Input parameters:
  - X,Y,W: Data points, function values, weight arrays (N points)
  - MAXDEG: Maximum fitting degree
  - EPS: Fitting precision control parameter
- Output parameters:
  - COEF: Polynomial coefficients (power series form)
  - NDEG: Actual fitting degree
  - R: Fitted values
- Key features:
  - Stable algorithm based on orthogonal polynomials
  - Supports weighted least squares fitting
  - Three degree selection modes: statistical testing, fixed degree, error control
  - F-test significance levels: 0.01, 0.05, 0.10
  - Automatic missing value handling
  - Numerically stable recursive relationships

## 130. smth9_dp.f [📝 src](fortran/smth9_dp.f)
**Function**: 9-point two-dimensional smoothing filter

**Detailed Description**:
- Contains one subroutine:
  - DSMTH9: 9-point weighted smoothing core algorithm
- Input parameters:
  - X: Input/output two-dimensional array (NI×NJ)
  - WRK: Work array (NI×NJ)
  - P,Q: Smoothing weight parameters (recommended P=0.5, Q=0.25)
  - XMSG: Missing value code
  - LWRAP: Boundary wraparound flag
- Smoothing formula:
  - $f_0' = f_0 + \frac{P}{4}\left(f_2 + f_4 + f_6 + f_8 - 4f_0\right) + \frac{Q}{4}\left(f_1 + f_3 + f_5 + f_7 - 4f_0\right)$
  - 9-point template: 8 neighboring points around center point
- Key features:
  - Based on Olson's modified 9-point smoothing algorithm
  - Supports periodic boundary conditions (longitude direction)
  - Automatically skips missing value points
  - Only smooths interior points, keeps boundaries unchanged
  - Suitable for meteorological field smoothing and noise filtering

## 131. sp2_util.f [📝 src](fortran/sp2_util.f)
**Function**: SPHEREPACK version 2 utility library (single precision)

**Detailed Description**:
- Contains multiple spherical coordinate transformation tools:
  - GEOMATV: Geographic coordinate vector field to mathematical coordinate transformation
  - MATGEOV: Mathematical coordinate vector field to geographic coordinate transformation
  - GEOMAT: Geographic coordinate scalar field to mathematical coordinate transformation
  - MATGEO: Mathematical coordinate scalar field to geographic coordinate transformation
  - TPLATX: Array transpose tool
  - CLATX: Latitude array flip tool
  - TRCTPR: Spherical harmonic coefficient truncation and tapering
  - GEOSCL: Array scalar multiplication
- Coordinate transformation characteristics:
  - Geographic coordinates: longitude×latitude, v component positive northward
  - Mathematical coordinates: latitude×longitude, v component negative
  - Includes array dimension transformation and latitude reversal
- Key features:
  - Based on John Adams' original code
  - Complete spherical coordinate system conversion
  - Supports vector and scalar field transformations
  - Includes spherical harmonic truncation and exponential tapering functions
  - Suitable for SPHEREPACK spherical harmonic analysis

## 132. sp2_util_dp.f [📝 src](fortran/sp2_util_dp.f)
**Function**: SPHEREPACK version 2 utility library (double precision version)

**Detailed Description**:
- Contains the same functionality as sp2_util.f, all in double precision versions:
  - DGEOMATV: Geographic to mathematical coordinate vector transformation
  - DMATGEOV: Mathematical to geographic coordinate vector transformation
  - DGEOMAT: Geographic to mathematical coordinate scalar transformation
  - DMATGEO: Mathematical to geographic coordinate scalar transformation
  - DTPLATX: Double precision array transpose
  - DCLATX: Double precision latitude reversal
  - DTRCTPR: Double precision spherical harmonic truncation
  - DGEOSCL: Double precision array scaling
- Transformation formulas:
  - Geographic→Mathematical: transpose + latitude reversal + v negation
  - Mathematical→Geographic: v negation + latitude reversal + transpose
- Key features:
  - Algorithm completely consistent with single precision version
  - Provides higher numerical precision
  - Suitable for high-precision spherical harmonic analysis calculations
  - Supports large-scale global model data processing

## 133. spareapoly.f [📝 src](fortran/spareapoly.f)
**Function**: Spherical polygon area calculation (Bevis-Cambareri algorithm)

**Detailed Description**:
- Contains three subroutines:
  - SPAREAPOLYI: Spherical polygon area calculation interface
  - SPAREAPOLY: Spherical polygon area core algorithm
  - TRNSFRMLON: Spherical coordinate transformation auxiliary function
- Input parameters:
  - VLAT,VLON: Polygon vertex latitudes and longitudes (degrees, NV vertices)
  - RAD: Sphere radius
- Output parameters:
  - AREA: Spherical polygon area (RAD$^2$ units)
- Theoretical basis:
  - Based on Bevis & Cambareri (1987) algorithm
  - Area formula: $A = \left(\sum\text{interior angles} - \pi(n-2)\right) \times R^2$
  - Vertices numbered in clockwise direction
- Key features:
  - Suitable for arbitrary shaped spherical polygons
  - Automatically handles duplicate vertices
  - Supports convex and concave polygons
  - Latitude $\pm$90°, longitude $\pm$180°
  - Suitable for GIS area calculations and geographic analysis

## 134. specx_dp.f [📝 src](fortran/specx_dp.f)
**Function**: Power spectral density estimation (double precision, complete statistical analysis)

**Detailed Description**:
- Contains three main subroutines:
  - DSPECX: Single variable power spectral analysis main interface
  - DSPECXY: Two-variable cross-spectral analysis main interface
  - DSPECXD: Spectral analysis core driver program
  - DTAPERX: Cosine bell tapering window function
  - DSWRNAV: Spectral estimation smoothing processing
- Input parameters:
  - X,Y: Time series data (NX points)
  - IOPT: Detrending options (0=remove mean, 1=linear, 2=quadratic)
  - JAVE: Number of spectral estimation smoothing points (odd number)
  - PCT: Tapering proportion [0-1]
  - SCL: Standardization options (-1,0,1,2 or user-defined)
- Output parameters:
  - FRQ,SPCX: Frequency and power spectrum
  - COSPC,QUSPC: Co-spectrum and quadrature spectrum
  - COHER,PHASE: Coherency and phase spectrum
  - SINFO: Statistical information (degrees of freedom, confidence intervals, etc.)
- Key features:
  - Based on FFT periodogram method
  - Modified Daniell smoothing algorithm
  - Complete statistical significance testing
  - Supports cosine bell tapering window
  - Suitable for climate time series frequency domain analysis

## 135. spi3.f [📝 src](fortran/spi3.f)
**Function**: Standardized Precipitation Index SPI-3 calculation (NCDC version)

**Detailed Description**:
- Contains multiple statistical analysis subroutines:
  - SPI3NCDC: NCL interface main program
  - SPIPE3: SPI core calculation algorithm
  - SAMLMR: L-moment parameter estimation
  - PELPE3: Pearson Type III parameter estimation
  - CDFPE3: Pearson Type III distribution function
  - QUASTN: Standard normal quantile function
- Input parameters:
  - PP: Monthly precipitation data (NTIM, starting from January)
  - NRUN: Rolling accumulation months (typical values 3,6,12,24)
  - PMSG: Missing value code
- Output parameters:
  - SPI3: Standardized precipitation index
  - PROBNE: Cumulative probability
  - PCPACC: Accumulated precipitation
- Algorithm characteristics:
  - Based on L-moment Pearson Type III distribution fitting
  - Mixed distribution handles zero precipitation frequency
  - SPI range constrained within $\pm$3.09
  - Calculate distribution parameters separately by month to maintain seasonality
- Key features:
  - US NCDC official SPI algorithm
  - Suitable for drought monitoring and climate analysis
  - Supports multiple time scales (3,6,12,24 months, etc.)
  - Includes complete statistical testing and L-moment estimation

## 136. spid.f [📝 src](fortran/spid.f)
**Function**: Standardized Precipitation Index SPI calculation (gamma distribution version)

**Detailed Description**:
- Contains multiple statistical analysis subroutines:
  - SPIGAMD: SPI gamma distribution main interface
  - ANVNRMD: Standard normal inverse transformation function
  - GAMFITD: Gamma distribution parameter estimation
  - GAMCDFD: Gamma distribution cumulative probability function
  - GAMINVD: Gamma distribution inverse function
  - Incomplete gamma function series: GAMMAPD, GAMMAQD, GAMSERD, GAMMCFD, etc.
- Input parameters:
  - PP: Precipitation time series data (NTIM)
  - NRUN: Accumulation time scale (months)
  - PMSG: Missing value code
- Output parameters:
  - INDEX: Standardized precipitation index values
- Mathematical principles:
  - Gamma distribution probability density: $f(x) = \frac{1}{\beta^\gamma \Gamma(\gamma)} x^{\gamma-1} e^{-x/\beta}$
  - Parameter estimation: $\alpha = \ln(\bar{x}) - \frac{\sum \ln(x)}{n}$, $\gamma = \frac{1 + \sqrt{1 + 4\alpha/3}}{4\alpha}$
  - Standardization transformation: $\text{SPI} = \Phi^{-1}[F_\gamma(x)]$
- Key features:
  - Based on Thom(1958) maximum likelihood estimation method
  - Supports mixed gamma distribution with zero values
  - Uses Abramowitz & Stegun normal inverse transformation
  - Fits distribution parameters separately by month
  - Suitable for drought monitoring and precipitation anomaly analysis

## 137. stattrop_dp.f [📝 src](fortran/stattrop_dp.f)
**Function**: WMO tropopause height calculation (double precision)

**Detailed Description**:
- Contains two subroutines:
  - STATTROPX: NCL interface wrapper program
  - DSTATTROP: WMO tropopause height core algorithm
- Input parameters:
  - PFULL,TFULL: Layer pressure and temperature (NLEV layers)
  - LAPSEC: Temperature lapse rate threshold (typical value 2 K/km)
  - PUNIT: Pressure units (0=hPa, 1=Pa)
- Output parameters:
  - PTROP: Tropopause pressure
  - ITROP: Tropopause layer index
- WMO definition criteria:
  - First tropopause: lowest layer with temperature lapse rate $\leq$ 2 K/km
  - 2 km average test: average lapse rate $\leq$ 2 K/km within 2 km above this layer
- Calculation method:
  - Temperature lapse rate: $\gamma = -\frac{dT}{dz} = \frac{d\ln T}{d\ln p} \times \frac{g}{R} \times 1000$
  - Hydrostatic approximation: $dz = -\frac{dp}{g\rho} = -\frac{RT}{gp}dp$
  - Logarithmic interpolation: $p_{\text{trop}} = \exp[p_1 + w(p_2 - p_1)]$
- Key features:
  - Based on WMO(1992) International Meteorological Vocabulary definition
  - Supports tropopause identification in 450-85 hPa range
  - Linear logarithmic pressure interpolation improves accuracy
  - Compatible with NCEP reanalysis tropopause algorithm
  - Suitable for atmospheric vertical structure analysis

## 138. statx.f [📝 src](fortran/statx.f)
**Function**: Statistical analysis toolkit (single precision, supports missing values)

**Detailed Description**:
- Contains complete statistical analysis subroutine set:
  - STAT2/STAT4: Second/fourth-order moment estimation
  - STAT2T: Trimmed mean and variance estimation
  - MEDMRNG: Median, range, and midrange estimation
  - ESAUTO/ESCROS: Autocorrelation and cross-correlation functions
  - VCVMNS/CORMNS: Variance-covariance and correlation matrices
  - ERRCHK: Outlier detection
  - SORTU/ISORTU: Quick sort algorithms
- Statistical quantity calculations:
  - Sample variance: $s^2 = \frac{\sum(x_i - \bar{x})^2}{n-1}$
  - Skewness coefficient: $\text{Skew} = \frac{\sum(x_i - \bar{x})^3}{ns^3}$
  - Kurtosis coefficient: $\text{Kurt} = \frac{\sum(x_i - \bar{x})^4}{ns^4} - 3$
  - Autocorrelation function: $r(\tau) = \frac{\sum(x_t - \bar{x})(x_{t+\tau} - \bar{x})}{\sum(x_t - \bar{x})^2}$
- Key features:
  - Comprehensive missing value handling mechanism
  - Trimmed statistics reduce extreme value influence
  - Numerically stable recursive algorithms
  - Symmetric storage mode saves memory
  - Includes complete matrix statistical operations
  - Suitable for time series and spatial data analysis

## 139. statx_dp.f [📝 src](fortran/statx_dp.f)
**Function**: Statistical analysis toolkit (double precision version)

**Detailed Description**:
- Contains double precision implementations of the same functionality as statx.f:
  - DSTAT2/DSTAT4: Double precision moment estimation
  - DSTAT2T: Double precision trimmed statistics
  - DMEDMRNG: Double precision median calculation
  - DESAUTO/DESCROS: Double precision correlation analysis
  - DVCVMNS/DCORMNS: Double precision matrix statistics
  - DSORTU/DISORTU: Double precision sorting algorithms
  - COLLAPSEXY: Missing value paired data extraction
- Numerical precision advantages:
  - Double precision floating point operations (approximately 15 significant digits)
  - Reduces accumulated rounding errors
  - Improves numerical stability for large sample statistics
  - Supports high-precision scientific computing requirements
- Key features:
  - Algorithm completely consistent with single precision version
  - Maintains same interface specifications
  - Suitable for high-precision statistical analysis
  - Numerical reliability for large dataset processing
  - Preferred choice for scientific computing and engineering applications

## 140. stndAtmUS76_dp.f [📝 src](fortran/stndAtmUS76_dp.f)
**Function**: US Standard Atmosphere 1976 thermodynamic calculations

**Detailed Description**:
- Contains three subroutines:
  - DSTDATMZ: Height → temperature, density, pressure
  - DSTDATMP: Pressure → height, temperature, density
  - STDZ2P/STDP2Z: Standard atmosphere core algorithms
- Standard atmospheric layer structure (0-84.852km):
  - Troposphere (0-11km): $T = 15.0 - 6.5z$ (°C)
  - Lower stratosphere (11-20km): $T = -56.5$ (°C)
  - Middle stratosphere (20-32km): $T = -56.5 + 1.0(z-20)$ (°C)
  - Upper stratosphere (32-47km): $T = -44.5 + 2.8(z-32)$ (°C)
- Thermodynamic relationships:
  - Hydrostatic equation: $dp = -\rho g dz$
  - Ideal gas: $p = \rho R T$, $R = 287.04$ J/(kg·K)
  - Isothermal layer: $p = p_0 \exp\left(-\frac{gz}{RT}\right)$
  - Linear gradient layer: $p = p_0 \left(\frac{T}{T_0}\right)^{-\frac{g}{R\gamma}}$
- Key features:
  - Based on US Standard Atmosphere 1976
  - Valid altitude range 0-84.852km
  - Pressure range 1013.25-0.00373hPa
  - Supports meter and feet unit conversion
  - High-precision double precision calculations
  - Suitable for aerospace and atmospheric physics applications

## 141. student_tprob.f [📝 src](fortran/student_tprob.f)
**Function**: Student's t-distribution probability calculation (based on incomplete beta function)

**Detailed Description**:
- Contains one subroutine:
  - STUPROBT: Student's t-distribution two-tailed probability calculation
- Input parameters:
  - T: t-statistic value
  - DF: Degrees of freedom (>0)
  - TMSG: Missing value code
- Output parameters:
  - RESULT: Two-tailed probability value P(|T|>t)
- Mathematical principles:
  - t-distribution probability density: $f(t) = \frac{\Gamma(\frac{\nu+1}{2})}{\sqrt{\nu\pi}\Gamma(\frac{\nu}{2})} \left(1+\frac{t^2}{\nu}\right)^{-\frac{\nu+1}{2}}$
  - Cumulative probability transformation: $P(|T|>t) = I_{\frac{\nu}{\nu+t^2}}\left(\frac{\nu}{2}, \frac{1}{2}\right)$
  - where $I_x(a,b)$ is the incomplete beta function
- Key features:
  - Based on SLATEC library incomplete beta function DBETAI
  - Avoids direct calculation of t-distribution integral through transformation
  - Automatically handles missing values and boundary conditions
  - Suitable for hypothesis testing and confidence interval calculations

## 142. svd_lap_dp.f [📝 src](fortran/svd_lap_dp.f)
**Function**: LAPACK singular value decomposition (double precision, based on Bretherton method)

**Detailed Description**:
- Contains three main subroutines:
  - DSVDLAP: SVD main interface (complete statistical analysis)
  - DSVDSV: SVD simplified interface (returns only singular vectors)
  - DSVDLAPACK1/DSVDLAPSV: SVD core algorithms
- Input parameters:
  - X,Y: Input data matrices (LDXY×NCX, LDXY×NCY)
  - MRT: Time dimension length
  - NCX,NCY: Spatial dimensions of X and Y
  - NSV: Number of singular values needed
  - IFLAG: Standardization option (1=standardize)
- Output parameters:
  - HOMLFT/HOMRGT: Left/right homogeneous correlation maps
  - HETLFT/HETRGT: Left/right heterogeneous correlation maps
  - AK,BK: Time expansion coefficients
  - PCVAR: Variance contribution percentage
- SVD decomposition process:
  - Cross-covariance matrix: $C = \frac{1}{T}X^T Y$
  - SVD decomposition: $C = U \Sigma V^T$
  - Time expansion coefficients: $A_k(t) = X(t) \cdot U_k$, $B_k(t) = Y(t) \cdot V_k$
- Key features:
  - Based on LAPACK DGESVD high-precision SVD algorithm
  - Complete Bretherton et al. statistical diagnostics
  - Supports data standardization and missing value handling
  - Includes condition number and rank estimation
  - Suitable for coupled mode analysis and dimensionality reduction

## 143. svd_util_dp.f [📝 src](fortran/svd_util_dp.f)
**Function**: SVD analysis utility toolkit (double precision)

**Detailed Description**:
- Contains four utility subroutines:
  - DSVDMTR: Matrix transpose tool
  - DSVDINFO: SVD information statistical analysis
  - DSVDPAR: Matrix data print output
  - DSVDPVC: Vector data print output
- DSVDINFO functionality:
  - Frobenius norm calculation: $||A||_F = \sqrt{\sum_{i=1}^{\min(m,n)} \sigma_i^2}$
  - Condition number estimation: $\kappa = \frac{\sigma_{\max}}{\sigma_{\min}}$
  - Numerical rank determination: $\text{rank}(A) = \#\{\sigma_i : \sigma_i > \text{TOL}\}$
  - Variance contribution: $\text{PC}_i = \frac{\sigma_i^2}{\sum \sigma_j^2} \times 100\%$
- Numerical stability:
  - Tolerance setting based on machine precision
  - Automatic singular value sorting
  - Robust rank estimation algorithm
- Key features:
  - Complete auxiliary tools for SVD analysis
  - Standardized output format for debugging
  - Double precision numerical accuracy guarantee
  - Modular design for easy maintenance

## 144. taper.f [📝 src](fortran/taper.f)
**Function**: Piecewise cosine bell tapering window function

**Detailed Description**:
- Contains one subroutine:
  - DTAPER: Time series tapering processing
- Input parameters:
  - X: Input time series (N points)
  - P: Tapering proportion [0,1] (e.g., P=0.1 means 10%)
  - IOPT: Tapering option (0=taper to series mean, 1=force taper to 0)
- Output parameters:
  - XT: Tapered time series
- Tapering window function:
  - Number of tapering points: $M = \max(1, \lfloor P \times N/2 + 0.5 \rfloor)$
  - Weight function: $w(i) = 0.5 - 0.5\cos\left(\frac{\pi(i-0.5)}{M}\right)$, $i = 1,2,\ldots,M$
  - Tapering formula: $x_t'(i) = (x(i) - \bar{x}) \times w(i) + \bar{x}$
- Application positions:
  - Series head: $i = 1, 2, \ldots, M$
  - Series tail: $i = N-M+1, N-M+2, \ldots, N$
  - Middle part remains unchanged
- Key features:
  - Based on Bloomfield Fourier analysis theory
  - Piecewise cosine bell smooth transition
  - Automatically detects constant series and forces tapering to 0
  - Reduces spectral leakage in FFT analysis
  - Suitable for frequency domain analysis of non-periodic signals

## 145. tdez1d.f [📝 src](fortran/tdez1d.f)
**Function**: NCAR Graphics three-dimensional trajectory marker plotting

**Detailed Description**:
- Contains one plotting subroutine:
  - TDEZ1D: Three-dimensional space trajectory marker visualization
- Input parameters:
  - X,Y,Z: Three-dimensional trajectory coordinate arrays (NX points)
  - IMRK: Marker type (1-5: tetrahedron, octahedron, cube, icosahedron, sphere)
  - RMRK: Marker radius
  - SMRK: Marker spacing (negative values allow overlap)
  - RMULT,THETA,PHI: Viewpoint spherical coordinates (multiplier, azimuth, elevation)
  - IST: Shading style index (1-8 color schemes)
- Viewpoint calculation:
  - Default position: $(R, \theta, \phi) = (2.5 \times DL, -55°, 70°)$
  - Cartesian coordinate transformation: 
    - $x_{eye} = x_{mid} + R \sin\phi \cos\theta$
    - $y_{eye} = y_{mid} + R \sin\phi \sin\theta$ 
    - $z_{eye} = z_{mid} + R \cos\phi$
  - Where $DL = \sqrt{(\Delta x)^2 + (\Delta y)^2 + (\Delta z)^2}$ is the diagonal length
- Color schemes:
  - IST=1: Wireframe mode
  - IST=2-8: Bottom grayscale + top color (red, green, blue, cyan, magenta, yellow)
  - Negative IST: White background + black foreground
- Key features:
  - Based on NCAR Graphics Tdpack three-dimensional plotting library
  - Supports multiple geometric markers and shading schemes
  - Automatically calculates optimal viewing angle and marker size
  - Triangulation rendering for realistic display
  - Suitable for three-dimensional trajectory and scatter data visualization

## 146. thornthwaite_v2.f [📝 src](fortran/thornthwaite_v2.f)
**Function**: Thornthwaite potential evapotranspiration calculation (version 2, supports missing values)

**Detailed Description**:
- Contains one subroutine:
  - THORN2: Thornthwaite potential evapotranspiration calculation main interface
- Input parameters:
  - T: Monthly mean temperature array (NTIM, °C)
  - LAT: Latitude (degrees)
  - TMSG: Missing value code
- Output parameters:
  - ETP: Potential evapotranspiration (NTIM, mm/month)
- Thornthwaite formula:
  - Heat index: $J = \sum_{i=1}^{12} \left(\frac{T_i}{5}\right)^{1.514}$ (only when $T_i > 0$)
  - Exponent coefficient: $c = 6.75 \times 10^{-7} J^3 - 7.71 \times 10^{-5} J^2 + 1.792 \times 10^{-2} J + 0.49239$
  - Evapotranspiration: $\text{ETP} = K \times 16 \times \left(\frac{10T}{J}\right)^c$ (when $T > 0$)
- Daylight correction coefficient:
  - Solar hour angle: $\omega = \arccos(-\tan\phi \tan\delta)$
  - Daylight hours: $N = \frac{24\omega}{\pi}$
  - Correction coefficient: $K = \frac{N}{12} \times \frac{\text{days}}{30}$
- Key features:
  - Based on original C code http://sac.csic.es/spei/spei_index.html
  - Uses monthly mean solar declination angle precomputed table
  - Automatically handles missing values and negative temperatures
  - Considers latitude and seasonal variations in daylight correction
  - Suitable for climatological evapotranspiration estimation and drought index calculations

## 147. triple2grid2d.f [📝 src](fortran/triple2grid2d.f)
**Function**: Nearest neighbor interpolation from scattered data to regular grid

**Detailed Description**:
- Contains one subroutine:
  - TRIPLE2GRID2D: Triple data to two-dimensional grid interpolation
- Input parameters:
  - X,Y,Z: Observation point coordinates and values (KPTS observation points)
  - LAT2D,LON2D: Target grid latitude and longitude (MDIM×NDIM)
  - DISTMX: Maximum search distance
  - MOPT: Distance calculation option (0=Cartesian coordinates, 1=great circle distance)
- Output parameters:
  - ZGRID: Interpolated grid data (MDIM×NDIM)
- Distance calculation formulas:
  - Cartesian coordinates: $d = \sqrt{(x_i - x_g)^2 + (y_i - y_g)^2}$
  - Great circle distance: $d = R_e \arccos[\sin\phi_g \sin\phi_i + \cos\phi_g \cos\phi_i \cos(\lambda_i - \lambda_g)]$
  - Where $R_e = 6371.22$ km is the Earth's radius
- Key features:
  - Nearest neighbor interpolation algorithm: each grid point is assigned the value of the nearest observation
  - Supports both Cartesian and spherical coordinate distance calculations
  - Distance threshold control: observation points beyond DISTMX do not participate in interpolation
  - Automatically handles missing value propagation
  - Suitable for fast interpolation from meteorological station data to grids

## 148. trssphx.f [📝 src](fortran/trssphx.f)
**Function**: Spherical harmonic scalar field grid transformation (SPHEREPACK interface)

**Detailed Description**:
- Contains one subroutine:
  - TRSSPHX: Scalar field spherical grid transformation main interface
- Input parameters:
  - INTL: Initialization flag (0=initialization, 1=transformation)
  - IGRIDA/IGRIDB: Source/target grid type ([$\pm$1,0/1]: equally-spaced/Gaussian, latitude/longitude priority)
  - NLONA,NLATA/NLONB,NLATB: Source/target grid dimensions
  - DA: Source grid scalar field data
  - JMODE: Truncation and tapering options
- Output parameters:
  - DB: Target grid scalar field data
- Spherical harmonic transformation process:
  - Analysis: $f(\lambda,\theta) \rightarrow \{A_n^m, B_n^m\}$
  - Truncation: Spherical harmonic coefficients truncated to specified wavenumber
  - Tapering: Optional exponential tapering to reduce Gibbs phenomenon
  - Synthesis: $\{A_n^m, B_n^m\} \rightarrow f'(\lambda',\theta')$
- Key features:
  - Based on NCAR SPHEREPACK double precision library
  - Supports equally-spaced and Gaussian grid interconversion
  - Supports different resolution grid transformations
  - Optional triangular truncation and exponential tapering
  - Maintains workspace-optimized efficient algorithms
  - Suitable for global model data post-processing

## 149. trvsphx.f [📝 src](fortran/trvsphx.f)
**Function**: Spherical harmonic vector field grid transformation (SPHEREPACK interface)

**Detailed Description**:
- Contains one subroutine:
  - TRVSPHX: Vector field spherical grid transformation main interface
- Input parameters:
  - INTL: Initialization flag (0=initialization, 1=transformation)
  - IGRIDA/IGRIDB: Source/target grid type
  - UA,VA: Source grid vector field components
  - IVECA/IVECB: Vector field type (0=mathematical coordinates, 1=geographical coordinates)
  - JMODE: Truncation and tapering options
- Output parameters:
  - UB,VB: Target grid vector field components
- Vector spherical harmonic transformation:
  - Analysis: $(u,v) \rightarrow \{A_n^m, B_n^m, C_n^m, D_n^m\}$
  - Where $(A,B)$ are vorticity potentials, $(C,D)$ are divergence potentials
  - Truncation and tapering processing applied to all coefficients
  - Synthesis: $\{A_n^m, B_n^m, C_n^m, D_n^m\} \rightarrow (u',v')$
- Coordinate system transformation:
  - Geographical coordinates: u eastward, v northward
  - Mathematical coordinates: u eastward, v southward (negative direction)
  - Automatically handles latitude ordering (north→south or south→north)
- Key features:
  - Based on NCAR SPHEREPACK vector spherical harmonic library
  - Preserves vector field divergence and vorticity properties
  - Supports coordinate system and grid type transformations
  - Optional truncation and tapering optimization
  - Suitable for wind field, ocean current and other vector field data processing

## 150. ttest.f [📝 src](fortran/ttest.f)
**Function**: Statistical hypothesis testing toolkit (t-test, F-test, correlation test)

**Detailed Description**:
- Contains four statistical testing subroutines:
  - DTTEST: Two-sample t-test (equal/unequal variances)
  - DFTEST: Two-sample F-test (variance equality)
  - DRTEST: Correlation coefficient significance test
  - DEQVSIZ: Effective sample size estimation (red noise process)
- t-test formulas:
  - Equal variances: $t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{s_p^2(1/n_1 + 1/n_2)}}$, $s_p^2 = \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2}$
  - Unequal variances (Welch): $t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{s_1^2/n_1 + s_2^2/n_2}}$, $\nu = \frac{(s_1^2/n_1 + s_2^2/n_2)^2}{(s_1^2/n_1)^2/(n_1-1) + (s_2^2/n_2)^2/(n_2-1)}$
- F-test formula:
  - $F = \frac{s_{\max}^2}{s_{\min}^2}$, degrees of freedom $\nu_1 = n_{\max} - 1$, $\nu_2 = n_{\min} - 1$
- Correlation test:
  - $t = r\sqrt{\frac{n-2}{1-r^2}}$, degrees of freedom $\nu = n-2$
  - Effective sample size: $n_{eff} = n \frac{1-r_1}{1+r_1}$ (when red noise is significant)
- Key features:
  - Based on SLATEC library incomplete beta function for p-value calculation
  - Complete hypothesis testing statistics and significance levels
  - Supports effective degrees of freedom correction for red noise time series
  - Covers most commonly used tests in climate statistics
  - Suitable for statistical significance analysis of meteorological data

## 151-201. 其余文件(超简版)  
## 151. upsl_dp.f [📝 src](fortran/upsl_dp.f)
**Function**: Sea level pressure calculation using hypsometric equation

**Detailed Description**:
- Contains two subroutines:
  - DPSLHY1: Sea level pressure calculation using virtual temperature
  - DPSLHY2: Sea level pressure calculation using temperature and mixing ratio
- Input parameters:
  - PRES: Lowest vertical level pressure (Pa)
  - Z: Lowest vertical level geopotential height (m)
  - TV: Virtual temperature (K) or T: Temperature (K)
  - W: Mixing ratio (kg/kg) (DPSLHY2)
  - XMSG: Missing value code
- Output parameters:
  - DPSLHY1/DPSLHY2: Sea level pressure (Pa)
- Key features:
  - Based on hypsometric equation: $p_{sl} = p \cdot \exp\left(\frac{gz}{R \cdot T_v}\right)$
  - Virtual temperature calculation in DPSLHY2: $T_v = T \cdot (1 + 0.608w)$
  - Uses standard gravitational acceleration: $g = 9.80665$ m/s$^2$
  - Dry air gas constant: $R = 287.04$ J/(kg·K)
  - Automatically handles missing values
  - Suitable for sea level pressure correction of meteorological station observation data

## 152. varimax_JiangLing_dp.f [📝 src](fortran/varimax_JiangLing_dp.f)
**Function**: Improved Varimax rotation algorithm for EOF analysis factor rotation

**Detailed Description**:
- Contains two subroutines:
  - ZROTEOF: Main rotation interface program, handles missing value detection and variance calculation
  - ZVORS: Core Varimax rotation algorithm (improved convergence criteria)
- Input parameters:
  - V: Eigenvector matrix (ND×NF)
  - NV,NF,ND: Number of variables, factors, array dimensions
  - EVAL: Eigenvalue array (NF)
  - VMSG: Missing value code
  - IOPT: Scaling option (0=no scaling, 1=scaling)
- Output parameters:
  - V: Rotated eigenvector matrix
  - ROTPCV: Rotated variance contribution percentage (NF)
  - ROTVAR: Rotated variance (NF)
  - KFLAG: Missing value presence flag
- Key features:
  - Based on Kaiser row normalization Varimax criterion
  - Improved convergence criteria: $|\gamma(k+1) - \gamma(k)| < 10^{-5}$
  - Rotation criterion: $\gamma = \sum_{j=1}^p \left[\sum_{i=1}^n a_{ij}^4 - \frac{1}{n}\left(\sum_{i=1}^n a_{ij}^2\right)^2\right]$
  - Supports robust handling of missing values
  - Suitable for identifying simple structure factors in climate mode analysis

## 153. varimax_dp.f [📝 src](fortran/varimax_dp.f)
**Function**: Classic Varimax factor rotation algorithm for principal component analysis

**Detailed Description**:
- Contains four subroutines:
  - ROTEOF: Main rotation interface program
  - VORS: Core Varimax rotation algorithm
  - DSUMFR: Vector sum/sum of squares calculation
  - DSCPFR: Vector scalar product calculation
  - VORSMSG: Missing value handling version rotation
- Input parameters:
  - V: Matrix to be rotated (ND×NF)
  - EVAL: Eigenvalues (NF)
  - IOPT: Processing option ($\leq$0=no scaling, >0=scaling)
- Output parameters:
  - V: Rotated factor loading matrix
  - ROTPCV: Rotated variance contribution percentage
  - ROTVAR: Rotated variance
- Key features:
  - Based on Veldman (1967) classic Varimax algorithm
  - Kaiser row normalization: $v_{ij}' = \frac{v_{ij}}{\sqrt{\sum_k v_{ik}^2}}$
  - Rotation angle: $\tan(4\theta) = \frac{D - 2AB/n}{C - (A^2-B^2)/n}$
  - Iterative convergence accuracy: $\epsilon = 0.0001$ radians
  - Compatible with IMSL FROTA function (w=1.0, norm=1, maxit=30)
  - Suitable for obtaining simple structure in factor analysis

## 154. vibeta_dp.f [📝 src](fortran/vibeta_dp.f)
**Function**: Vertical integration using Kuo/Trenberth diagnostic method

**Detailed Description**:
- Contains four subroutines:
  - DVIBETA: Main vertical integration algorithm
  - DINT2P2: Pressure coordinate interpolation
  - DINITVC: Vector initialization
  - DCOPYVC: Vector copying
- Input parameters:
  - P: Pressure level array (NLEV, Pa)
  - X: Array to be integrated (NLEV)
  - PSFC: Surface pressure (Pa)
  - PBOT,PTOP: Integration lower and upper boundary pressures (Pa)
  - LINLOG: Interpolation type (1=linear, 2=logarithmic)
  - XMSG: Missing value code
- Output parameters:
  - VINT: Vertical integration result
  - IER: Error code
- Key features:
  - Based on Kuo & Trenberth (1991) diagnostic method
  - Integration formula: $\int_{P_{top}}^{P_{bot}} X \frac{dp}{g} = \sum_{k} \beta_k X_k \Delta p_k$
  - Weight coefficient: $\beta_k = \frac{p_{sfc} - p_{k+1}}{p_{k-1} - p_{k+1}}$ (surface layer)
  - Supports sounding data processing with missing values
  - Suitable for atmospheric column total calculations (water vapor, ozone, etc.)

## 155. vinth2p_dp.f [📝 src](fortran/vinth2p_dp.f)
**Function**: Vertical interpolation from CCM2/3 hybrid coordinate data to pressure coordinates

**Detailed Description**:
- Contains two subroutines:
  - VINTH2P: Main vertical interpolation program
  - PRNT: Debug printing tool
- Input parameters:
  - DATI: Hybrid coordinate data (IMAX×NLAT×NLEVI)
  - HBCOFA,HBCOFB: Hybrid coordinate coefficients A and B (NLEVIP1)
  - P0: Reference pressure (mb)
  - PLEVO: Output pressure levels (NLEVO, mb)
  - PSFC: Surface pressure (IMAX×NLAT, Pa)
  - INTYP: Interpolation type (1=linear, 2=logarithmic, 3=double logarithmic)
  - KXTRP: Extrapolation flag (0=use special values, 1=extrapolate)
- Output parameters:
  - DATO: Pressure coordinate data (IMAX×NLAT×NLEVO)
- Key features:
  - Hybrid coordinate formula: $P(k) = A(k) \times P_0 + B(k) \times P_{sfc}$
  - Linear interpolation: $f = f_1 + (f_2-f_1) \frac{p-p_1}{p_2-p_1}$
  - Logarithmic interpolation: $f = f_1 + (f_2-f_1) \frac{\ln(p/p_1)}{\ln(p_2/p_1)}$
  - Double logarithmic interpolation: $f = f_1 + (f_2-f_1) \frac{A2LN(p)-A2LN(p_1)}{A2LN(p_2)-A2LN(p_1)}$
  - Where $A2LN(x) = \ln(\ln(x+2.72))$
  - Supports boundary extrapolation and special value handling
  - Suitable for atmospheric model vertical coordinate transformation

## 156. vinth2p_ecmwf.f [📝 src](fortran/vinth2p_ecmwf.f)
**Function**: Vertical interpolation from ECMWF model hybrid coordinate data to pressure coordinates (with ECMWF extrapolation)

**Detailed Description**:
- Contains one subroutine:
  - VINTH2PECMWF: Main vertical interpolation program (with ECMWF extrapolation algorithm)
- Input parameters:
  - DATI: Hybrid coordinate data (IMAX×NLAT×NLEVI)
  - HBCOFA,HBCOFB: Hybrid coordinate coefficients A and B (NLEVIP1)
  - P0: Reference pressure (mb)
  - PLEVO: Output pressure levels (NLEVO, mb)
  - PSFC: Surface pressure (IMAX×NLAT, Pa)
  - VARFLG: Variable flag (-1=geopotential height, +1=temperature, 0=other)
  - TBOT: Lowest level temperature (IMAX×NLAT, K)
  - PHIS: Surface geopotential (IMAX×NLAT, m$^2$/s$^2$)
  - INTYP: Interpolation type (1=linear, 2=logarithmic, 3=double logarithmic)
- Output parameters:
  - DATO: Pressure coordinate data (IMAX×NLAT×NLEVO)
- Key features:
  - Based on ECMWF IFS model extrapolation algorithm
  - Temperature extrapolation: $T = T^* \left(1 + \alpha \ln p + \frac{1}{2}\alpha^2(\ln p)^2 + \frac{1}{6}\alpha^3(\ln p)^3\right)$
  - Geopotential height extrapolation: $Z = Z_s - \frac{RT^*}{g}\ln\left(\frac{p}{p_s}\right)\left(1 + \frac{1}{2}\alpha\ln p + \frac{1}{6}\alpha^2(\ln p)^2\right)$
  - Where $\alpha = \frac{0.0065 \times R}{g}$ is the standard atmosphere lapse rate parameter
  - Supports complex height correction algorithm (2000-2500m transition layer)
  - Suitable for high-precision vertical interpolation of ECMWF data

## 157. vinth2p_ecmwf_nodes_dp.f [📝 src](fortran/vinth2p_ecmwf_nodes_dp.f)
**Function**: Vertical interpolation of ECMWF hybrid coordinate node data (one-dimensional node array version)

**Detailed Description**:
- Contains one subroutine:
  - DVINTH2PECMWFNODES: Node-version ECMWF vertical interpolation algorithm
- Input parameters:
  - DATI: Hybrid coordinate data (NPTS×NLEVI)
  - HBCOFA,HBCOFB: Hybrid coordinate coefficients A and B
  - PSFC: Surface pressure (NPTS, Pa)
  - TBOT: Lowest level temperature (NPTS, K)
  - PHIS: Surface geopotential (NPTS, m$^2$/s$^2$)
  - VARFLG: Variable type flag
- Output parameters:
  - DATO: Interpolated data (NPTS×NLEVO)
- Key features:
  - Suitable for efficient processing of one-dimensional node arrays
  - Uses the same ECMWF extrapolation algorithm
  - Temperature correction formula: $T^* = T_{bot} \times \left(1 + \alpha\left(\frac{p_{sfc}}{p_{lev}} - 1\right)\right)$
  - Height-related temperature gradient processing: standard gradient below 2000m, transition gradient 2000-2500m
  - Double precision calculation ensures numerical accuracy
  - Suitable for observation points or irregular grid data interpolation

## 158. vinth2p_nodes_dp.f [📝 src](fortran/vinth2p_nodes_dp.f)
**Function**: Standard hybrid coordinate node data vertical interpolation (without ECMWF extrapolation)

**Detailed Description**:
- Contains one subroutine:
  - DVINTH2PNODES: Standard node-version vertical interpolation algorithm
- Input parameters:
  - DATI: Hybrid coordinate data (NPTS×NLEVI)
  - HBCOFA,HBCOFB: Hybrid coordinate coefficients
  - PLEVO: Output pressure levels (NLEVO, mb)
  - PSFC: Surface pressure (NPTS, Pa)
  - KXTRP: Extrapolation flag (0=use special values, 1=simple extrapolation)
- Output parameters:
  - DATO: Interpolated data (NPTS×NLEVO)
- Key features:
  - Does not contain ECMWF complex extrapolation algorithm
  - Three interpolation methods: linear, logarithmic, double logarithmic
  - Double logarithmic function: $A2LN(x) = \ln(\ln(x + 2.72))$
  - Simple extrapolation: uses lowest level data values
  - Hybrid coordinate formula: $P(k) = A(k) \times P_0 + B(k) \times P_{sfc}$
  - Suitable for standard atmospheric model data interpolation

## 159. vintp2p_ecmwf.f [📝 src](fortran/vintp2p_ecmwf.f)
**Function**: ECMWF model pressure coordinate data to pressure coordinate vertical interpolation

**Detailed Description**:
- Contains one subroutine:
  - VINTP2PECMWF: Pressure to pressure interpolation (with ECMWF extrapolation)
- Input parameters:
  - DATI: Input pressure coordinate data (IMAX×NLAT×NLEVI)
  - PRESI: Input pressure field (IMAX×NLAT×NLEVI, Pa)
  - PLEVO: Output pressure levels (NLEVO, mb)
  - PSFC: Surface pressure (IMAX×NLAT, Pa)
  - VARFLG: Variable type (-1=height, +1=temperature, 0=other)
  - TBOT,PHIS: Surface temperature and geopotential
- Output parameters:
  - DATO: Interpolated pressure coordinate data (IMAX×NLAT×NLEVO)
- Key features:
  - Based on interpolation algorithm using existing pressure field
  - Uses the same ECMWF extrapolation formula as vinth2p_ecmwf
  - Three interpolation types: $f = f_1 + (f_2-f_1) \times \frac{I(p)-I(p_1)}{I(p_2)-I(p_1)}$
  - Where $I(p)$ is the interpolation function: linear, $\ln p$ or $\ln(\ln(p+2.72))$
  - Suitable for reanalysis data and observational data vertical interpolation
  - Does not require hybrid coordinate coefficients, directly uses pressure field

## 160. volAve.f [📝 src](fortran/volAve.f)
**Function**: Calculate weighted volume average of three-dimensional data fields

**Detailed Description**:
- Contains two subroutines:
  - DWGTVOLAVE: Volume averaging using separated weights
  - DWGTVOLAVECCM: Volume averaging using 3D weights (suitable for CCM models)
- Input parameters:
  - T: 3D data field (MX×NY×KZ)
  - WGTX: X-direction weights (MX)
  - WGTY: Y-direction weights (NY)  
  - WGTZ: Z-direction weights (KZ) or 3D weights (MX×NY×KZ)
  - IFLAG: Missing value handling flag (0=ignore, 1=return missing value when encountered)
- Output parameters:
  - AVE: Volume weighted average
- Key features:
  - Two-step weighted averaging algorithm:
    1. Vertical weighting: $\bar{T}_z(x,y) = \frac{\sum_k T(x,y,k) \cdot w_z(k)}{\sum_k w_z(k)}$
    2. Horizontal weighting: $\langle T \rangle = \frac{\sum_{x,y} \bar{T}_z(x,y) \cdot w_x(x) \cdot w_y(y)}{\sum_{x,y} w_x(x) \cdot w_y(y)}$
  - CCM version supports 3D weights (such as $\Delta p$ field): $w_z(x,y,k) = \Delta p(x,y,k)$
  - Automatically handles missing values
  - Suitable for atmospheric column mass-weighted averaging, global averaging calculations  
## 161. volRmse.f [📝 src](fortran/volRmse.f)
**Function**: Calculate weighted volume root mean square error between two three-dimensional data fields

**Detailed Description**:
- Contains two subroutines:
  - DWGTVOLRMSE: Three-dimensional RMSE calculation using separated weights
  - DWGTVOLRMSECCM: CCM version RMSE calculation using 3D weights
- Input parameters:
  - T,Q: Two 3D data fields (MX×NY×KZ)
  - WGTX: X-direction weights (MX)
  - WGTY: Y-direction weights (NY) 
  - WGTZ: Z-direction weights (KZ) or 3D weights (MX×NY×KZ)
  - DPT,DPQ: 3D weight fields in CCM version (such as layer thickness Δp)
  - IFLAG: Missing value handling flag
- Output parameters:
  - RMSE: Volume weighted root mean square error
- Key features:
  - Two-step RMSE calculation:
    1. Vertical RMSE: $\text{RMSE}_z(x,y) = \sqrt{\frac{\sum_k w_z(k) \cdot (T(x,y,k)-Q(x,y,k))^2}{\sum_k w_z(k)}}$
    2. Horizontal weighting: $\text{RMSE} = \sqrt{\frac{\sum_{x,y} w_x(x) \cdot w_y(y) \cdot \text{RMSE}_z^2(x,y)}{\sum_{x,y} w_x(x) \cdot w_y(y)}}$
  - CCM version uses average layer thickness: $w_{avg}(x,y,k) = \frac{1}{2}(w_T(x,y,k) + w_Q(x,y,k))$
  - Suitable for three-dimensional model validation, data comparison, inter-model error analysis

## 162. vpower.f [📝 src](fortran/vpower.f)
**Function**: Calculate space-time power spectra, cross-spectra, coherency and phase of two variables

**Detailed Description**:
- Contains four subroutines:
  - SPCTIMCROSS1: Main control program, calculates array dimensions and calls core algorithms
  - SPCTIMCROSS2: Core space-time spectral analysis algorithm
  - SPCTIMCROSS3: Smoothing processing and coherency calculation
  - SMTH121STCP: 1-2-1 smoothing filter
- Input parameters:
  - X,Y: Two 3D input fields (NL×NM×NT), latitude×longitude×time
  - NL,NM,NT: Number of latitude points, longitude points, time points
- Output parameters:
  - STC: Space-time spectral results (NLP1×NT2P1×16), containing 16 spectral components
- Algorithm flow:
  - Spatial FFT: $\tilde{X}(n,t) = \text{FFT}_{\text{lon}}[X(\text{lon},\text{lat},t)]$
  - Temporal FFT: $\hat{X}(n,\omega) = \text{FFT}_{\text{time}}[\tilde{X}(n,t)]$
  - Power spectrum: $P_{XX}(n,\omega) = |\hat{X}(n,\omega)|^2$
  - Cross-spectrum: $P_{XY}(n,\omega) = \text{Re}[\hat{X}^*(n,\omega)\hat{Y}(n,\omega)]$
  - Coherency: $\text{Coh}^2(n,\omega) = \frac{P_{XY}^2 + Q_{XY}^2}{P_{XX} \cdot P_{YY}}$
- Key features:
  - Distinguishes symmetric and antisymmetric components (relative to equator)
  - 1-2-1 smoothing filter reduces noise
  - Calculates phase relationships and vector components
  - Suitable for atmospheric wave analysis, Wheeler-Kiladis diagram generation

## 163. wavelet.f [📝 src](fortran/wavelet.f)
**Function**: Continuous wavelet transform analysis, supporting multiple mother wavelets and significance testing

**Detailed Description**:
- Contains multiple subroutines:
  - WAVELETI: NCL interface main program
  - WAVELET: Core wavelet transform algorithm
  - WAVE_FUNCTION: Mother wavelet function calculation
  - WAVE_SIGNIF: Significance testing
  - STATSWAVE: Time series statistics calculation
- Supported mother wavelets:
  - Morlet wavelet (mother=0): $\psi(t) = \pi^{-1/4}e^{ik_0t}e^{-t^2/2}$
  - Paul wavelet (mother=1): $\psi(t) = \frac{2^m i^m m!}{\sqrt{\pi(2m)!}}(1-it)^{-(m+1)}$  
  - DOG wavelet (mother=2): $\psi(t) = \frac{(-1)^{m+1}}{\sqrt{\Gamma(m+1/2)}}D^m[e^{-t^2/2}]$
- Input parameters:
  - Y: Time series data (N points)
  - DT: Time interval
  - MOTHER: Mother wavelet type (0-2)
  - PARAM: Wavelet parameters (k₀, m, etc.)
  - S0,DJ,JTOT: Scale parameters
- Output parameters:
  - WAVE: Wavelet coefficients (N×JTOT×2)
  - POWER: Wavelet power spectrum
  - PHASE: Wavelet phase spectrum  
  - SIGNIF: Significance levels
- Key features:
  - Based on Torrence & Compo (1998) algorithm
  - Complete statistical significance testing
  - Cone of influence (COI) calculation
  - Global wavelet spectrum (GWS)
  - Suitable for time series periodic analysis, climate oscillation research

## 164. wetbulb_profs.f [📝 src](fortran/wetbulb_profs.f)
**Function**: Calculate wet-bulb temperature profiles based on Stipanuk (1973) algorithm

**Detailed Description**:
- Contains multiple subroutines:
  - WETBULBPROFS: Main interface program, handles profile data
  - TWBPROFS: Core wet-bulb temperature calculation
  - OXPROFS,WXPROFS: Potential temperature and mixing ratio calculation
  - ESATPROFS: Saturated water vapor pressure calculation
  - TSAPROFS,TDAPROFS: Moist adiabatic and dry adiabatic line temperatures
- Input parameters:
  - T: Temperature profile (°C, NPTS points)
  - TD: Dew point temperature profile (°C, NPTS points)
  - P: Pressure profile (mb, NPTS points)
- Output parameters:
  - TWB: Wet-bulb temperature profile (°C, NPTS points)
- Calculation steps:
  - Mixing ratio line: $w = \frac{622 \cdot e_s(T_d)}{p - e_s(T_d)}$
  - Dry adiabatic line: $\theta = T \cdot \left(\frac{1000}{p}\right)^{0.286}$
  - Iterative intersection pressure: $p_i = p \cdot 2^{0.02(T_{mr}(w,p_i) - T_{da}(\theta,p_i))}$
  - Moist adiabatic line: $\theta_e = T \cdot \left(\frac{1000}{p}\right)^{0.286} \cdot \exp\left(\frac{2.65 \cdot w}{T}\right)$
- Key features:
  - Based on NOAA PROFS system algorithm
  - Uses Nordquist (1973) saturated water vapor pressure formula
  - Automatically handles missing values
  - High-precision iterative algorithm (10 iterations)
  - Suitable for atmospheric physics calculations, weather forecasting, radiosonde data processing

## 165. wgtVertAvg_beta.f [📝 src](fortran/wgtVertAvg_beta.f)
**Function**: Mass-weighted vertical integration or averaging using Beta factors

**Detailed Description**:
- Contains three subroutines:
  - DWVBETAP1: One-dimensional pressure coordinate version
  - DWVBETAP3: Three-dimensional pressure coordinate version
  - DWVBETAP: Core Beta factor calculation algorithm
- Input parameters:
  - P: Pressure level array (KLEV or MLON×NLAT×KLEV)
  - X: Variable to be integrated (MLON×NLAT×KLEV)
  - PSFC: Surface pressure (MLON×NLAT)
  - IPUNIT: Pressure unit flag (0=mb, 1=Pa)
  - IOPT: Output option (0=integration, 1=average)
  - PTOP,PBOT: Integration upper and lower bounds
- Output parameters:
  - XVB: Vertical integration or averaging results (MLON×NLAT)
- Beta factor calculation:
  - Layer thickness: $\Delta p(k) = p(k+1) - p(k-1)$
  - Beta factor: $\beta(k) = \frac{\min(p_{bot}, p_{sfc}) - p(k-1)}{p(k+1) - p(k-1)}$
  - Integration: $I = \sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)$
  - Average: $\bar{X} = \frac{\sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)}{\sum_k \beta(k) \cdot \Delta p(k)}$
- Key features:
  - Based on Trenberth (1991) mass-conserving diagnostic method
  - Automatically handles boundary layer and terrain effects
  - Supports partial layer integration (PTOP to PBOT)
  - Suitable for atmospheric column total calculations, mass-weighted averaging, vertical integration diagnostics  
## 166. wk_smooth121.f [📝 src](fortran/wk_smooth121.f)
**Function**: 1-2-1 weighted smoothing filter

**Detailed Description**:
- Contains one subroutine:
  - WKSMOOTH121: 1-2-1 weight smoothing filter core algorithm
- Input parameters:
  - VV: Input/output data array (VN)
  - VN: Array length
  - NN: Number of valid data points ($\leq$VN)
  - SPV: Missing value identifier
  - DUM: Work array (VN)
- Output parameters:
  - VV: Smoothed data array (in-place modification)
- Key features:
  - Smoothing weight formula: $f'(i) = \frac{1 \cdot f(i-1) + 2 \cdot f(i) + 1 \cdot f(i+1)}{4}$
  - Boundary handling (first point): $f'(1) = \frac{3 \cdot f(1) + 1 \cdot f(2)}{4}$
  - Boundary handling (last point): $f'(n) = \frac{1 \cdot f(n-1) + 3 \cdot f(n)}{4}$
  - Conservative sum property: $\sum f'(i) = \sum f(i)$
  - Automatically skips missing values
  - Suitable for time series smoothing in Wheeler-Kiladis diagram analysis

## 167. wk_utils.f [📝 src](fortran/wk_utils.f)
**Function**: Wheeler-Kiladis diagram analysis utility toolkit

**Detailed Description**:
- Contains three subroutines:
  - WKSMOOTH121: 1-2-1 smoothing filter (same as wk_smooth121.f)
  - WKTAPERTOZERO: Time series tapering to zero processing
  - WKDETREND: Linear trend removal
- WKTAPERTOZERO input parameters:
  - TS: Time series data (N)
  - N,NMI,NN,TP: Array length, starting index, valid length, taper points
- WKDETREND input parameters:
  - X2: Time series data (NX, in-place modification)
  - NX,VMI,NV: Array length, starting position, valid data length
- Key features:
  - Tapering window function: $w(j) = 0.5 \times (1 - \cos(\frac{(j-1)\pi}{tp}))$, $j = 1,2,\ldots,tp$
  - Tapering formula: $ts'(i) = ts(i) \times w(i)$
  - Linear detrending: $X'(t) = X(t) - (a + b \cdot t)$
  - Where $b = \frac{\sum t \cdot X(t) - n\bar{t}\bar{X}}{\sum t^2 - n\bar{t}^2}$, $a = \bar{X} - b\bar{t}$
  - Tapering processing satisfying FFT periodicity requirements
  - Suitable for atmospheric wave spectral analysis preprocessing

## 168. wrf_bint3d.f [📝 src](fortran/wrf_bint3d.f)
**Function**: WRF model coordinate transformation and three-dimensional bilinear interpolation

**Detailed Description**:
- Contains multiple subroutines:
  - DMAPTFORM: WRF map projection coordinate transformation
  - DBINT3D: Three-dimensional bilinear interpolation main program
  - DBINT: Two-dimensional bilinear interpolation algorithm
  - DONED: One-dimensional interpolation core function
- DMAPTFORM input parameters:
  - DSKMC,MIYCORS,MJXCORS: Grid spacing and dimensions
  - NPROJ,XLATC,XLONC: Projection type, central latitude and longitude
  - TRUE1,TRUE2: True latitude parameters
  - RIY,RJX,RLAT,RLON: Grid coordinates and geographic coordinates
  - IDIR: Conversion direction (1=grid→geographic, -1=geographic→grid)
- DBINT3D input parameters:
  - DATA_IN: Input three-dimensional data (NX×NY×NZ)
  - OBSII,OBSJJ: Observation point grid coordinates (NOBSICRS×NOBSJCRS)
  - ICRS,JCRS: Cross flag parameters
- Output parameters:
  - DATA_OUT: Interpolated data (NOBSICRS×NOBSJCRS×NZ)
- Key features:
  - Supports three map projections:
    - Lambert Conformal (NPROJ=1): $\text{cone} = \frac{\ln(\cos(\phi_1)/\cos(\phi_2))}{\ln(\tan(\pi/4-\phi_2/2)/\tan(\pi/4-\phi_1/2))}$
    - Polar Stereographic (NPROJ=2): $\text{cone} = 1$
    - Mercator (NPROJ=3): $y = R_e \ln(\tan(\pi/4 + \phi/2))$
  - Four-point bilinear interpolation: $f(x,y) = \sum_{i,j} w_{ij} f_{ij}$
  - Weight functions: $w_{ij} = (1-u)(1-v), u(1-v), (1-u)v, uv$
  - Automatically handles grid boundaries and missing values
  - Suitable for WRF model data post-processing and observation point interpolation

## 169. wrf_cloud_fracf.f90 [📝 src](fortran/wrf_cloud_fracf.f90)
**Function**: WRF model cloud fraction calculation (Fortran90 implementation)

**Detailed Description**:
- Contains two subroutines:
  - DCLOUDFRAC: Fixed pressure threshold cloud fraction calculation (deprecated)
  - DCLOUDFRAC2: Flexible vertical coordinate cloud fraction calculation (recommended)
- DCLOUDFRAC2 input parameters:
  - RH: Relative humidity field (EW×NS×NZ, %)
  - VERT: Vertical coordinate field (EW×NS×NZ, pressure or height)
  - LOW_THRESH,MID_THRESH,HIGH_THRESH: Low, middle, high cloud thresholds
  - VERT_INC_W_HEIGHT: Vertical coordinate direction flag (0=decreasing, 1=increasing)
  - MSG: Missing value code
- Output parameters:
  - LOWC,MIDC,HIGHC: Low, middle, high cloud fractions (EW×NS)
- Key features:
  - Cloud layer classification standards:
    - Low clouds: $p > 970$ hPa or $z < 2$ km
    - Middle clouds: $700 < p < 970$ hPa or $2 < z < 6$ km  
    - High clouds: $p < 450$ hPa or $z > 10$ km
  - Cloud fraction calculations:
    - Low clouds: $CF_{low} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - Middle clouds: $CF_{mid} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - High clouds: $CF_{high} = \max(0, \min(1, 2.5 \times \frac{\text{RH}_{max}}{100} - 1.5))$
  - Uses OpenMP parallelization for accelerated computation
  - Supports arbitrary vertical coordinate systems
  - Automatically handles terrain masking situations
  - Suitable for WRF model cloud diagnostics and climate analysis

## 170. wrf_constants.f90 [📝 src](fortran/wrf_constants.f90)
**Function**: WRF model physical constants definition module

**Detailed Description**:
- Contains one Fortran90 module:
  - WRF_CONSTANTS: Complete physical constants definition collection
- Main physical constants:
  - Basic constants:
    - $\pi = 3.141592653589793$
    - Earth radius: $R_e = 6.37 \times 10^6$ m
    - Gravitational acceleration: $g = 9.81$ m/s$^2$
  - Thermodynamic constants:
    - Dry air gas constant: $R_d = 287$ J/(kg·K)
    - Water vapor gas constant: $R_v = 461.6$ J/(kg·K)
    - Specific heat at constant pressure: $C_p = 1004.5$ J/(kg·K)
    - Adiabatic exponent: $\gamma = R_d/C_p = 0.286$
  - Water vapor parameters:
    - Molecular weight ratio: $\epsilon = \frac{M_d}{M_v} = 0.622$
    - Saturation vapor pressure constant: $e_0 = 6.112$ hPa
    - Clausius-Clapeyron constants: $L_1 = 17.67$, $L_2 = 29.65$
  - Cloud physics parameters:
    - Cloud water absorption coefficient: $\alpha_w = 0.145$ m$^2$/g
    - Cloud ice absorption coefficient: $\alpha_i = 0.272$ m$^2$/g
    - Water density: $\rho_w = 1000$ kg/m$^3$
- Key features:
  - Fully compatible with WRF model module_model_constants.F
  - Double precision floating point definitions ensure numerical accuracy
  - Contains complete missing value definition system
  - Supports default fill values for multiple data types
  - Moist adiabatic process correction parameters:
    - $C_{p,moist} = C_p(1 + 0.887 \times q_v)$
    - $R_{gas,moist} = R_d(1 + 0.608 \times q_v)$
  - Suitable for standardized constant references in all WRF-related computation programs  
## 171. wgtVertAvg_beta.f [📝 src](fortran/wgtVertAvg_beta.f)
**Function**: Mass-weighted vertical integration or averaging using Beta factors

**Detailed Description**:
- Contains three subroutines:
  - DWVBETAP1: One-dimensional pressure coordinate version
  - DWVBETAP3: Three-dimensional pressure coordinate version  
  - DWVBETAP: Core Beta factor calculation algorithm
- Input parameters:
  - P: Pressure level array (KLEV or MLON×NLAT×KLEV)
  - X: Variable to be integrated (MLON×NLAT×KLEV)
  - PSFC: Surface pressure (MLON×NLAT)
  - IPUNIT: Pressure unit flag (0=mb, 1=Pa)
  - IOPT: Output option (0=integration, 1=average)
  - PTOP,PBOT: Integration upper and lower bounds
- Output parameters:
  - XVB: Vertical integration or averaging results (MLON×NLAT)
- Key features:
  - Beta factor calculation: $\beta(k) = \frac{\min(p_{bot}, p_{sfc}) - p(k-1)}{p(k+1) - p(k-1)}$
  - Layer thickness: $\Delta p(k) = p(k+1) - p(k-1)$
  - Integration formula: $I = \sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)$
  - Average formula: $\bar{X} = \frac{\sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)}{\sum_k \beta(k) \cdot \Delta p(k)}$
  - Based on Trenberth (1991) mass-conserving diagnostic method
  - Automatically handles boundary layer and terrain effects
  - Suitable for atmospheric column total calculations and mass-weighted averaging

## 172. wk_smooth121.f [📝 src](fortran/wk_smooth121.f)
**Function**: 1-2-1 weighted smoothing filter

**Detailed Description**:
- Contains one subroutine:
  - WKSMOOTH121: 1-2-1 weight smoothing filter core algorithm
- Input parameters:
  - VV: Input/output data array (VN)
  - VN: Array length
  - NN: Number of valid data points ($\leq$VN)  
  - SPV: Missing value identifier
  - DUM: Work array (VN)
- Output parameters:
  - VV: Smoothed data array (in-place modification)
- Key features:
  - Smoothing weight formula: $f'(i) = \frac{1 \cdot f(i-1) + 2 \cdot f(i) + 1 \cdot f(i+1)}{4}$
  - Boundary handling (first point): $f'(1) = \frac{3 \cdot f(1) + 1 \cdot f(2)}{4}$
  - Boundary handling (last point): $f'(n) = \frac{1 \cdot f(n-1) + 3 \cdot f(n)}{4}$
  - Conservative sum property: $\sum f'(i) = \sum f(i)$
  - Automatically skips missing values
  - Suitable for time series smoothing in Wheeler-Kiladis diagram analysis

## 173. wk_utils.f [📝 src](fortran/wk_utils.f)
**Function**: Wheeler-Kiladis diagram analysis utility toolkit

**Detailed Description**:
- Contains three subroutines:
  - WKSMOOTH121: 1-2-1 smoothing filter
  - WKTAPERTOZERO: Time series tapering to zero processing
  - WKDETREND: Linear trend removal
- WKTAPERTOZERO input parameters:
  - TS: Time series data (N)
  - N,NMI,NN,TP: Array length, starting index, valid length, taper points
- WKDETREND input parameters:
  - X2: Time series data (NX, in-place modification)
  - NX,VMI,NV: Array length, starting position, valid data length
- Key features:
  - Tapering window function: $w(j) = 0.5 \times (1 - \cos(\frac{(j-1)\pi}{tp}))$, $j = 1,2,\ldots,tp$
  - Tapering formula: $ts'(i) = ts(i) \times w(i)$
  - Linear detrending: $X'(t) = X(t) - (a + b \cdot t)$
  - Slope: $b = \frac{\sum t \cdot X(t) - n\bar{t}\bar{X}}{\sum t^2 - n\bar{t}^2}$, $a = \bar{X} - b\bar{t}$
  - Tapering processing satisfying FFT periodicity requirements
  - Suitable for atmospheric wave spectral analysis preprocessing


## 176. wrf_constants.f90 [📝 src](fortran/wrf_constants.f90)
**Function**: WRF model physical constants definition module

**Detailed Description**:
- Contains one Fortran90 module:
  - WRF_CONSTANTS: Complete physical constants definition collection
- Main physical constants:
  - Basic constants:
    - $\pi = 3.1415926535897932$
    - Earth radius: $R_e = 6.37 \times 10^6$ m
    - Gravitational acceleration: $g = 9.81$ m/s$^2$
  - Thermodynamic constants:
    - Dry air gas constant: $R_d = 287$ J/(kg·K)
    - Water vapor gas constant: $R_v = 461.6$ J/(kg·K)
    - Specific heat at constant pressure: $C_p = 1004.5$ J/(kg·K)
    - Adiabatic exponent: $\gamma = R_d/C_p = 0.286$
  - Water vapor parameters:
    - Molecular weight ratio: $\epsilon = \frac{M_d}{M_v} = 0.622$
    - Saturation vapor pressure constant: $e_0 = 6.112$ hPa
    - Clausius-Clapeyron constants: $L_1 = 17.67$, $L_2 = 29.65$
  - Cloud physics parameters:
    - Cloud water absorption coefficient: $\alpha_w = 0.145$ m$^2$/g
    - Cloud ice absorption coefficient: $\alpha_i = 0.272$ m$^2$/g
    - Water density: $\rho_w = 1000$ kg/m$^3$
- Key features:
  - Fully compatible with WRF model module_model_constants.F
  - Double precision floating point definitions ensure numerical accuracy
  - Contains complete missing value definition system
  - Supports default fill values for multiple data types
  - Moist adiabatic process correction parameters:
    - $C_{p,moist} = C_p(1 + 0.887 \times q_v)$
    - $R_{gas,moist} = R_d(1 + 0.608 \times q_v)$
  - Suitable for standardized constant references in all WRF-related computation programs

## 177. wrf_fctt.f90 [📝 src](fortran/wrf_fctt.f90)
**Function**: WRF cloud top temperature calculation (based on optical thickness method)

**Detailed Description**:
- Contains one subroutine:
  - WRFCTTCALC: Cloud top temperature calculation main program
- Input parameters:
  - PRS: Pressure field (ew×ns×nz, Pa)
  - TK: Temperature field (ew×ns×nz, K)
  - QCI,QCW: Cloud ice and cloud water mixing ratios (ew×ns×nz, g/kg)
  - QVP: Water vapor mixing ratio (ew×ns×nz, g/kg)
  - GHT: Geopotential height (ew×ns×nz, m)
  - TER: Terrain height (ew×ns, m)
  - HAVEQCI: Whether cloud ice data exists (0/1)
  - FILL_NOCLOUD: No-cloud fill option (0/1)
  - OPT_THRESH: Optical thickness threshold (typical value 1.0)
- Output parameters:
  - CTT: Cloud top temperature (ew×ns, K)
  - PF: Full layer pressure field (ew×ns×nz, Pa)
- Key features:
  - Optical thickness calculations:
    - Cloud water: $\tau_w = \alpha_w \times q_{cw} \times \frac{\Delta p}{g}$
    - Cloud ice: $\tau_i = \alpha_i \times q_{ci} \times \frac{\Delta p}{g}$
    - Total optical thickness: $\tau = \tau_w + \tau_i$
  - Cloud top detection: Integrate downward from model top until $\tau \geq 1.0$
  - Temperature interpolation: Linear interpolation to obtain cloud top temperature
  - Uses OpenMP parallelization for accelerated computation
  - Automatically handles terrain and boundary conditions
  - Suitable for satellite data validation and cloud microphysics diagnostics

## 178. wrf_pvo.f90 [📝 src](fortran/wrf_pvo.f90)
**Function**: WRF potential vorticity and absolute vorticity calculation

**Detailed Description**:
- Contains two subroutines:
  - DCOMPUTEABSVORT: Absolute vorticity calculation
  - DCOMPUTEPV: Potential vorticity calculation
- DCOMPUTEABSVORT input parameters:
  - U,V: Wind field components (nxp1×ny×nz, nx×nyp1×nz, m/s)
  - MSFU,MSFV,MSFT: Map scale factors
  - COR: Coriolis parameter (nx×ny, s$^{-1}$)
  - DX,DY: Grid spacing (m)
- DCOMPUTEPV input parameters:
  - U,V: Wind field components
  - THETA: Potential temperature field (nx×ny×nz, K)
  - PRS: Pressure field (nx×ny×nz, Pa)
  - Map projection and grid parameters
- Output parameters:
  - AV: Absolute vorticity (nx×ny×nz, ×10⁵ s$^{-1}$)
  - PV: Potential vorticity (nx×ny×nz, ×10⁻$^2$ PVU)
- Key features:
  - Absolute vorticity formula: $\zeta_{abs} = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + f$
  - Potential vorticity formula: $PV = -g(\frac{\partial \theta}{\partial p}\zeta_{abs} - \frac{\partial v}{\partial p}\frac{\partial \theta}{\partial x} + \frac{\partial u}{\partial p}\frac{\partial \theta}{\partial y})$
  - Map projection correction: Uses map scale factor $m$ for coordinate transformation
  - Finite differences: Uses central differences to calculate partial derivatives
  - Boundary handling: Automatically handles one-sided differences at grid boundaries
  - Uses OpenMP parallelization
  - Suitable for dynamical diagnostics and weather analysis

## 179. wrf_pw.f90 [📝 src](fortran/wrf_pw.f90)  
**Function**: WRF precipitable water calculation

**Detailed Description**:
- Contains one subroutine:
  - DCOMPUTEPW: Atmospheric column precipitable water calculation
- Input parameters:
  - P: Pressure field (nx×ny×nz, Pa)
  - TV: Virtual temperature field (nx×ny×nz, K)
  - QV: Water vapor mixing ratio (nx×ny×nz, kg/kg)
  - HT: Geopotential height (nx×ny×nzh, m), nzh=nz+1
- Output parameters:
  - PW: Precipitable water (nx×ny, kg/m$^2$ or mm)
- Key features:
  - Integration formula: $PW = \sum_{k=1}^{nz} \frac{p(k)}{R_d \times T_v(k)} \times q_v(k) \times \Delta h(k)$
  - Where $\Delta h(k) = h(k+1) - h(k)$ is the layer thickness
  - Density calculation: $\rho = \frac{p}{R_d \times T_v}$
  - Mass integration: $PW = \sum \rho \times q_v \times \Delta h$
  - Unit conversion: kg/m$^2$ = mm water column height
  - Uses OpenMP parallelization
  - Vertical integration from surface to model top
  - Suitable for atmospheric water vapor analysis and precipitation forecasting

## 180. wrf_relhl.f90 [📝 src](fortran/wrf_relhl.f90)
**Function**: WRF storm-relative helicity calculation

**Detailed Description**:
- Contains one subroutine:
  - DCALRELHL: Storm-relative helicity calculation main program
- Input parameters:
  - U,V: Wind field components (miy×mjx×mkzh, m/s)
  - GHT: Geopotential height (miy×mjx×mkzh, m)
  - TER: Terrain height (miy×mjx, m)
  - LAT: Latitude (miy×mjx, degrees)
  - TOP: Integration upper limit height (m, typical value 3000m)
- Output parameters:
  - SREH: Storm-relative helicity (miy×mjx, m$^2$/s$^2$)
- Key features:
  - Storm motion speed estimation:
    - Mean wind speed: $\vec{U}_{avg} = \frac{\sum_{k} \vec{U}(k) \Delta h(k)}{\sum_{k} \Delta h(k)}$ (3-10km layer average)
    - Storm velocity: $\vec{C} = 0.75 |\vec{U}_{avg}|$, direction 30° to the right (Northern Hemisphere)
  - Helicity calculation: $SRH = -\sum_{k} [(\vec{U}(k) - \vec{C}) \times (\vec{U}(k) - \vec{U}(k+1))]_z$
  - Physical meaning: Measures streamline vorticity in convective storm inflow environment
  - Threshold references:
    - SRH < 100 m$^2$/s$^2$: Cutoff value
    - SRH = 150-299: Supercells possibly producing weak tornadoes
    - SRH = 300-499: Favorable for strong tornado development
    - SRH > 450: Violent tornadoes
  - Uses OpenMP parallelization
  - Automatically identifies Northern/Southern Hemispheres
  - Suitable for severe convective weather forecasting and tornado warning  

## 181. wrf_write_wps.f [📝 src](fortran/wrf_write_wps.f)
**Function**: WRF preprocessing system intermediate format file writing function

**Detailed Description**:
- Contains one subroutine:
  - WRITE_INTERMEDIATE_WPS: Write meteorological data to WPS intermediate format file
- Input parameters:
  - OUTPUT_NAME: Output file name prefix
  - FIELDIN: Variable name (length $\leq$9 characters)
  - UNITSIN: Variable units (length $\leq$25 characters)  
  - DESCIN: Variable description (length $\leq$46 characters)
  - HDATEIN: Date-time string (format YYYY-MM-DD_HH:MM:SS)
  - DATA_SOURCE: Data source identifier (length $\leq$32 characters)
  - XLVL: Vertical layer/pressure level value
  - IPROJ: Map projection type (0-5)
  - START_LOCATION: Starting position identifier (length $\leq$8 characters)
  - STARTLAT,STARTLON: Starting latitude and longitude
  - DELTALAT,DELTALON: Latitude and longitude intervals
  - XLONC: Projection center longitude  
  - TRUELAT1,TRUELAT2: Projection true latitudes
  - NLATS: Number of Gaussian projection latitudes
  - DX,DY: Grid spacing (m)
  - NX,NY: Grid dimensions
  - IS_WIND_EARTH_REL: Wind field relative to Earth coordinate flag
  - IFV: WPS format version number (default 5)
  - XFCST: Forecast time (hours)
  - EARTH_RADIUS: Earth radius (km)
  - SLAB: Two-dimensional data field (NX×NY)
- Supported map projection types:
  - IPROJ=0: Cylindrical equidistant projection (lat-lon grid)
  - IPROJ=1: Mercator projection
  - IPROJ=3: Lambert conformal conic projection  
  - IPROJ=4: Gaussian projection
  - IPROJ=5: Polar stereographic projection
- Output file format:
  - File name format: `$\text{OUTPUT\_NAME}:\text{YYYY-MM-DD\_HH:MM:SS}$`
  - Binary big-endian format
  - WPS intermediate format record structure:
    - Record 1: Format version number IFV
    - Record 2: General information (date, variable name, units, description, etc.)
    - Record 3: Projection-related parameters (varies with IPROJ)
    - Record 4: Wind field coordinate system flag
    - Record 5: Data field SLAB
- Key features:
  - Supports file append mode, multiple variables for the same time can be written to the same file
  - Automatically handles parameter writing for different projections
  - Compatible with WPS v3.0+ format standards
  - Uses big-endian to ensure cross-platform compatibility
  - Suitable for preprocessing data preparation for WRF numerical weather prediction models  
## 182. writematrix.f [📝 src](fortran/writematrix.f)
**Function**: Multi-data type matrix formatted file I/O function library

**Detailed Description**:
- Contains five subroutines:
  - WRITEMATRIXI: Integer matrix writing function
  - WRITEMATRIXF: Single precision real matrix writing function
  - WRITEMATRIXD: Double precision real matrix writing function
  - WRITEMATRIXB: Byte-type integer matrix writing function (INTEGER*1)
  - WRITEMATRIXS: Short integer matrix writing function (INTEGER*2)
- Input parameters:
  - FNAME: Output file name (if "*" then output to standard output)
  - NROW,NCOL: Number of matrix rows and columns
  - X: Input matrix data (NCOL×NROW)
  - FMTX: Fortran format descriptor string
  - TITLE: Output title string
  - TITSP: Number of spaces before title
  - IOPT: Output option flag
- Output format options:
  - IOPT=0: Output matrix data only
    - Format: $\text{FORMAT} = \text{"("} || \text{FMTX} || \text{")"}$
  - IOPT=1: Output row index + matrix data
    - Format: $\text{FORMAT} = \text{"(I5,1X,"} || \text{FMTX} || \text{")"}$
- Title formatting:
  - When TITSP>0: $\text{TITFMT} = \text{"("} || \text{TITSP} || \text{"X,A)"}$
  - When TITSP$\leq$0: $\text{TITFMT} = \text{"(X,A)"}$
- Supported data types:
  - INTEGER: Standard integer type
  - REAL: Single precision floating point
  - DOUBLE PRECISION: Double precision floating point  
  - INTEGER*1: Byte-type integer (-128 to 127)
  - INTEGER*2: Short integer (-32768 to 32767)
- Key features:
  - Supports both file output and standard output modes
  - Flexible format control, supports user-defined Fortran format specifiers
  - Optional row index output (numbered from 0)
  - Controllable title output position and format
  - Matrix output in row-major order: $\text{Row i} = [X(1,i), X(2,i), \ldots, X(\text{NCOL},i)]$
  - Automatic buffer flushing ensures output completeness
  - Supports non-standard Fortran extension data types (INTEGER*1/2)
  - Suitable for formatted output of scientific computing results and data exchange  

## 183. wrunave_dp.f [📝 src](fortran/wrunave_dp.f)
**Function**: Weighted running average filtering algorithm (double precision)

**Detailed Description**:
- Contains two subroutines:
  - DWGTRUNAVE: Weighted running average main interface
  - DWRUNAVX77: Weighted running average core algorithm implementation
- Input parameters:
  - X: Input time series data (NPTS, in-place modification)
  - NPTS: Number of data points
  - WGT: Weight vector (NWGT)
  - NWGT: Moving window size (number of weight points)
  - KOPT: Boundary handling option
  - XMSG: Missing value code
  - WORK: Work array (length $\geq$NPTS+2×NHALF, NHALF=NWGT/2)
  - LWORK: Work array length
- Output parameters:
  - X: Filtered time series (in-place modification)
  - IER: Error code (0=success, <0=error)
- Boundary handling options:
  - KOPT<0: Cyclic boundary conditions
    - Example (NWGT=3): $X(1) = w_1 X(n) + w_2 X(1) + w_3 X(2)$
    - Example (NWGT=4): $X(1) = w_1 X(n-1) + w_2 X(n) + w_3 X(1) + w_4 X(2)$
  - KOPT=0: Boundary points set to missing values
    - Points at both ends insufficient for window size set to XMSG
    - Example (NWGT=3): $X(1) = \text{XMSG}$, $X(2) = w_1 X(1) + w_2 X(2) + w_3 X(3)$
  - KOPT>0: Reflective (symmetric) boundary conditions  
    - Example (NWGT=3): $X(1) = w_1 X(2) + w_2 X(1) + w_3 X(2)$
    - Example (NWGT=4): $X(1) = w_1 X(3) + w_2 X(2) + w_3 X(1) + w_4 X(2)$
- Weighted average formula:
  - Weight normalization: $\text{WSUM} = \sum_{i=1}^{\text{NWGT}} w_i$
  - If WSUM>1: $\text{WSUM} = \frac{1}{\text{WSUM}}$ (normalization factor)
  - Filtering formula: $X'(n) = \text{WSUM} \times \sum_{k=0}^{\text{NWGT}-1} w_{k+1} \times X(n-\text{NHALF}+k)$
- Missing value handling:
  - When any point within window is missing, output is set to XMSG
  - Automatically skip missing value points, not participating in filtering calculation
- Error handling:
  - IER=-11: NPTS$\leq$0
  - IER=-12: NWGT>NPTS (window size exceeds data length)
- Key features:
  - Double precision floating point computation ensures numerical accuracy
  - Supports arbitrary length weight vectors
  - Three boundary handling strategies adapt to different application needs
  - Automatic weight normalization ensures filter gain
  - Memory-optimized work array design
  - In-place modification saves memory space
  - Suitable for time series smoothing, noise suppression, signal filtering  

## 184. xy1pdf77.f [📝 src](fortran/xy1pdf77.f)
**Function**: One-dimensional probability density function estimation algorithm (Fortran77 implementation)

**Detailed Description**:
- Contains one subroutine:
  - X1PDF77: Calculate probability density function (histogram) for one-dimensional data
- Input parameters:
  - NX: Total number of data points
  - X: Input data array (NX)
  - XMSG: Missing value code
  - NBX: Number of histogram bins
  - NBXP1: Number of bin boundary points (NBX+1)
  - BINXBND: Bin boundary array (NBXP1)
  - IPCNT: Output unit flag (0=frequency, 1=percentage)
- Output parameters:
  - PDF: Probability density function values (NBX)
  - IER: Error code (0=success, 1=no valid data)
- Algorithm flow:
  1. Initialize all PDF bins to 0
  2. Count valid data points: $K_x = \sum_{i=1}^{NX} \mathbf{1}_{X(i) \neq \text{XMSG}}$
  3. Data binning process:
     - For bin $j$ (j=1,2,...,NBX): 
     - Counting condition: $\text{BINXBND}(j) \leq X(i) < \text{BINXBND}(j+1)$
     - Special handling for last boundary: when $X(i) = \text{BINXBND}(\text{NBX}+1)$ assign to last bin
  4. Frequency calculation: $\text{PDF}(j) = \sum_{i=1}^{NX} \mathbf{1}_{X(i) \in \text{bin }j}$
- Output unit conversion:
  - IPCNT=0: Output raw frequencies
  - IPCNT=1: Output percentages $\text{PDF}(j) = 100 \times \frac{\text{PDF}(j)}{K_x}$
- Performance optimization features:
  - Missing value detection optimization: skip missing value check when $K_x = NX$
  - Double-loop structure adapted for large dataset processing
- Boundary handling:
  - Left boundary inclusive: $X(i) \geq \text{BINXBND}(j)$
  - Right boundary exclusive: $X(i) < \text{BINXBND}(j+1)$ 
  - Rightmost boundary special inclusion: $X(i) = \text{BINXBND}(\text{NBX}+1)$
- Key features:
  - Supports arbitrary user-defined bin boundaries
  - Automatic missing value handling
  - Can output both frequency and percentage formats
  - Performance optimization for large datasets
  - Standard histogram statistical algorithm
  - Suitable for data distribution analysis, statistical modeling, data visualization preprocessing  

## 185. xy2pdf77.f [📝 src](fortran/xy2pdf77.f)
**Function**: Two-dimensional joint probability density function estimation algorithm (Fortran77 implementation)

**Detailed Description**:
- Contains one subroutine:
  - XY2PDF77: Calculate joint probability density function for two-dimensional data (2D histogram)
- Input parameters:
  - NXY: Total number of data pairs
  - X,Y: Input data arrays (NXY)
  - XMSG,YMSG: Missing value codes for X and Y
  - NBY: Number of histogram bins in Y direction
  - MBX: Number of histogram bins in X direction
  - MBXP1: Number of bin boundary points in X direction (MBX+1)
  - NBYP1: Number of bin boundary points in Y direction (NBY+1)
  - BINXBND: X direction bin boundary array (MBXP1)
  - BINYBND: Y direction bin boundary array (NBYP1)
  - IPCNT: Output unit flag (0=frequency, 1=percentage)
- Output parameters:
  - PDF: Two-dimensional probability density function matrix (MBX×NBY)
  - IER: Error code (0=success, 1=no valid data pairs)
- Algorithm flow:
  1. Initialize all PDF bins to 0
  2. Count valid data pairs: $K_{xy} = \sum_{i=1}^{NXY} \mathbf{1}_{X(i) \neq \text{XMSG} \land Y(i) \neq \text{YMSG}}$
  3. Two-dimensional data binning process:
     - For bin $(m,n)$ (m=1,...,MBX; n=1,...,NBY):
     - X counting condition: $\text{BINXBND}(m) \leq X(i) < \text{BINXBND}(m+1)$
     - Y counting condition: $\text{BINYBND}(n) \leq Y(i) < \text{BINYBND}(n+1)$
  4. Joint frequency calculation: $\text{PDF}(m,n) = \sum_{i=1}^{NXY} \mathbf{1}_{(X(i),Y(i)) \in \text{bin }(m,n)}$
- Output unit conversion:
  - IPCNT=0: Output raw frequencies
  - IPCNT=1: Output percentages $\text{PDF}(m,n) = 100 \times \frac{\text{PDF}(m,n)}{K_{xy}}$
- Missing value handling strategy:
  - Strict pairing: data pair is valid only when $X(i) \neq \text{XMSG}$ and $Y(i) \neq \text{YMSG}$
  - If either variable is missing, the entire data pair is excluded
- Performance optimization features:
  - Missing value detection optimization: skip missing value check when $K_{xy} = NXY$
  - Triple nested loop structure: outer double loop iterates through bins, inner loop iterates through data
- Boundary handling:
  - Left boundary inclusive: $X(i) \geq \text{BINXBND}(m)$, $Y(i) \geq \text{BINYBND}(n)$
  - Right boundary exclusive: $X(i) < \text{BINXBND}(m+1)$, $Y(i) < \text{BINYBND}(n+1)$
- Key features:
  - Supports arbitrary user-defined two-dimensional bin grids
  - Strict missing value paired processing
  - Can output both frequency and percentage formats
  - Performance optimization for large datasets
  - Standard two-dimensional histogram statistical algorithm
  - Suitable for 2D data distribution analysis, correlation studies, joint probability modeling  

## 186. xy2pdf90.f [📝 src](fortran/xy2pdf90.f)
**Function**: Two-dimensional joint probability density function estimation algorithm (Fortran90 optimized implementation)

**Detailed Description**:
- Contains one subroutine:
  - XY2PDF90: Calculate joint probability density function for two-dimensional data (Fortran90 vectorized version)
- Input parameters:
  - NXY: Total number of data pairs
  - X,Y: Input data arrays (NXY)
  - XMSG,YMSG: Missing value codes for X and Y
  - NBY: Number of histogram bins in Y direction
  - MBX: Number of histogram bins in X direction
  - MBXP1: Number of bin boundary points in X direction (MBX+1)
  - NBYP1: Number of bin boundary points in Y direction (NBY+1)
  - BINXBND: X direction bin boundary array (MBXP1)
  - BINYBND: Y direction bin boundary array (NBYP1)
  - IPCNT: Output unit flag (0=frequency, 1=percentage)
- Output parameters:
  - PDF: Two-dimensional probability density function matrix (MBX×NBY)
  - IER: Error code (0=success, 1=no valid data pairs)
- Algorithm flow:
  1. Fortran90 array initialization: $\text{PDF} = 0.0$
  2. Vectorized valid data counting: $K_{xy} = \text{COUNT}(X \neq \text{XMSG} \land Y \neq \text{YMSG})$
  3. Vectorized two-dimensional binning process:
     - For bin $(m,n)$ (m=1,...,MBX; n=1,...,NBY):
     - Vectorized condition: $\text{LOGICAL MASK} = (X \geq \text{BINXBND}(m)) \land (X < \text{BINXBND}(m+1))$ $\land (Y \geq \text{BINYBND}(n)) \land (Y < \text{BINYBND}(n+1))$
  4. Vectorized frequency calculation: $\text{PDF}(m,n) = \text{COUNT}(\text{LOGICAL MASK})$
- Output unit conversion:
  - IPCNT=0: Output raw frequencies
  - IPCNT=1: Vectorized percentage conversion $\text{PDF} = 100 \times \frac{\text{PDF}}{K_{xy}}$
- Fortran90 optimization features:
  - Built-in COUNT function: replaces explicit loop counting, improves execution efficiency
  - Array assignment syntax: $\text{PDF} = 0.0$ and $\text{PDF} = \text{UNIT} \times (\text{PDF}/K_{xy})$
  - Logical array operations: vectorization using .AND. and .NE. logical operators
  - Avoids explicit inner loops: uses COUNT intrinsic function to replace triple nested loops
- Performance improvements:
  - Reduces number of loop levels compared to Fortran77 version
  - Utilizes compiler vectorization optimization
  - Optimized memory access patterns
- Missing value handling strategy:
  - Strict pairing: data pair is valid only when $X(i) \neq \text{XMSG}$ and $Y(i) \neq \text{YMSG}$
  - Vectorized missing value detection and filtering
- Boundary handling:
  - Left boundary inclusive: $X(i) \geq \text{BINXBND}(m)$, $Y(i) \geq \text{BINYBND}(n)$
  - Right boundary exclusive: $X(i) < \text{BINXBND}(m+1)$, $Y(i) < \text{BINYBND}(n+1)$
- Key features:
  - Fortran90 modern syntax and intrinsic functions
  - Vectorized operations significantly improve performance
  - Concise code structure and higher readability
  - Maintains completely consistent algorithm logic with Fortran77 version
  - Automatic compiler optimization support
  - Suitable for large-scale 2D data distribution analysis, high-performance statistical computing  

## 187. z2geouv_dp.f [📝 src](fortran/z2geouv_dp.f)
**Function**: Geopotential height to geostrophic wind conversion algorithm (double precision)

**Detailed Description**:
- Contains three subroutines:
  - Z2GEOUV: Main interface program, automatically determines latitude ordering
  - ZUVNEW: North-to-south data reordering processing program  
  - Z2GUV: Geostrophic wind core calculation program (requires south-to-north ordering)
- Input parameters:
  - Z: Geopotential height field (MLON×NLAT, units: gpm)
  - MLON,NLAT: Number of longitude and latitude grid points
  - ZMSG: Missing value code
  - GLON: Longitude array (MLON)
  - GLAT: Latitude array (NLAT)
  - IOPT: Boundary condition option (0=non-periodic, 1=periodic)
- Output parameters:
  - UG,VG: Geostrophic wind U and V components (MLON×NLAT, units: m/s)
- Basic geostrophic wind equations:
  - Geostrophic balance: $f \vec{V_g} = -\nabla \Phi \times \hat{k}$
  - U component: $U_g = -\frac{g}{f} \frac{\partial Z}{\partial y} = -\frac{g}{f} \frac{1}{R} \frac{\partial Z}{\partial \phi}$
  - V component: $V_g = \frac{g}{f} \frac{\partial Z}{\partial x} = \frac{g}{f} \frac{1}{R \cos\phi} \frac{\partial Z}{\partial \lambda}$
- Physical constants and parameters:
  - Gravitational acceleration: $g = 9.80616 \text{ m/s}^2$
  - Earth radius: $R_e = 6371220 \text{ m}$
  - Earth rotation angular velocity: $\Omega = 7.292 \times 10^{-5} \text{ rad/s}$
  - Coriolis parameter: $f = 2\Omega \sin\phi$
- Calculation methods:
  - Gravity/Coriolis force ratio: $\frac{g}{f} = \frac{g}{2\Omega \sin\phi}$
  - Spatial difference distances: 
    - Zonal: $\Delta x = 2R \Delta\lambda \cos\phi$ (units: m)
    - Meridional: $\Delta y = R \Delta\phi$ (units: m)
  - Central difference scheme: $\frac{\partial Z}{\partial x} \approx \frac{Z(i+1,j) - Z(i-1,j)}{2\Delta x}$
- Boundary handling strategies:
  - Latitude boundaries: use one-sided differences, distance weights halved
  - Longitude boundaries: 
    - IOPT=1: Periodic boundary conditions $Z(0,j) = Z(\text{MLON},j)$
    - IOPT=0: Copy adjacent point values, multiply V component by 2 to compensate for one-sided difference
- Special region handling:
  - Near equator: when $|\phi| < 10^{-5}$, Coriolis term set to missing value
  - Near poles: when $|\phi| > 89.999°$, zonal distance set to 0
- Data ordering requirements:
  - Algorithm requires latitude ordering from south to north
  - Automatically detects input data ordering
  - North-to-south data automatically reordered
- Key features:
  - Double precision floating point calculation ensures accuracy
  - Automatic data ordering detection and processing
  - Flexible boundary condition settings
  - Complete set of geophysical parameters
  - Accurate spherical geometry calculations
  - Suitable for numerical weather prediction, atmospheric dynamics analysis, wind field diagnostics

## 188. zon_mpsi.f [📝 src](fortran/zon_mpsi.f)
**Function**: Zonal mean stream function calculation algorithm (atmospheric circulation analysis)

**Detailed Description**:
- Contains two subroutines:
  - DZPSIDRV: Main driver program, allocates work arrays
  - DZONMPSI: Zonal mean stream function core calculation program
- Input parameters:
  - V: Meridional wind field (MLON×NLAT×KLEV, units: m/s)
  - LAT: Latitude array (NLAT, units: degrees)
  - P: Pressure level array (KLEV, units: Pa)
  - PS: Surface pressure (MLON×NLAT, units: Pa)  
  - VMSG: Missing value code
  - MLON,NLAT,KLEV: Number of longitude, latitude, and vertical levels
- Output parameters:
  - ZMPSI: Zonal mean stream function (NLAT×KLEV, units: kg·m/s$^2$)
- Basic stream function calculation equations:
  - Continuity equation: $\frac{\partial \bar{v}}{\partial p} + \frac{\partial \bar{\omega}}{\partial y} = 0$
  - Stream function definition: $\bar{v} = -\frac{\partial \Psi}{\partial p}$, $\bar{\omega} = \frac{1}{a\cos\phi}\frac{\partial \Psi}{\partial \phi}$
  - Integral relationship: $\Psi(\phi,p) = \int_{p_{top}}^{p} \bar{v}(\phi,p') dp'$
- Calculation methods and parameters:
  - Geophysical constants:
    - Gravitational acceleration: $g = 9.80616 \text{ m/s}^2$
    - Earth radius: $a = 6.37122 \times 10^6 \text{ m}$
    - Conversion constant: $C = \frac{2\pi a}{g}$
  - Latitude weight: $C(\phi) = C \times \cos\phi$
- Algorithm flow:
  1. Pressure level extension: create half-layer pressures $P_{half}(k) = \frac{P(k)+P(k+1)}{2}$
  2. Pressure difference calculation: $\Delta P(k) = P(k+1) - P(k-1)$
  3. Meridional wind field preprocessing: set grid points with $P > P_s$ to missing values
  4. Zonal mean calculation: $\bar{v}(\phi,p) = \frac{1}{N_{lon}} \sum_{i=1}^{N_{lon}} v(i,\phi,p)$
  5. Vertical integration: $\Psi(\phi,p) = -C(\phi) \sum_{k=1}^{p} \bar{v}(\phi,k) \Delta P(k)$
- Boundary condition handling:
  - Top boundary: $P_{top} = 500 \text{ Pa}$ (lower stratosphere boundary)
  - Bottom boundary: $P_{bot} = 100500 \text{ Pa}$ (below surface)
  - Lower boundary condition: $\Psi(\phi,P_{bot}) = 0$ (mass conservation requirement)
- Physical meaning of stream function:
  - Positive values: clockwise circulation (Northern Hemisphere)
  - Negative values: counterclockwise circulation (Northern Hemisphere)
  - Contour lines: represent streamlines of mean meridional circulation
- Numerical processing characteristics:
  - Automatic missing value propagation: missing values at any level affect integration results
  - Mass conservation correction: force bottom boundary stream function to zero
  - CSM convention: output values take negative sign to conform to climate model convention
- Key features:
  - Rigorous derivation based on mass continuity equation
  - Supports arbitrary vertical coordinate systems
  - Automatically handles regions below topography
  - Complete atmospheric circulation diagnostics
  - Suitable for Hadley circulation and Walker circulation analysis
  - Supports standardized processing of climate model output
  - Suitable for atmospheric circulation diagnostics, climate change research, meridional heat transport analysis

## 189. zregr.f [📝 src](fortran/zregr.f)
**Function**: Multi-variable linear regression analysis algorithm (with standardized coefficients)

**Detailed Description**:
- Contains two subroutines:
  - DZREGR1: Main multi-variable regression program, calculates standardized coefficients
  - DZREGR2: Core regression algorithm, least squares solution
- Input parameters:
  - Y: Dependent variable time series (N)
  - X: Independent variable matrix (N×M)
  - N: Number of observation samples
  - M: Number of independent variables
  - XMSG,YMSG: Missing value codes for X and Y
  - Work arrays: WK,YY,COV,XSD,XMEAN,A,AINV,S
- Output parameters:
  - C: Original regression coefficients (M)
  - CNORM: Standardized regression coefficients (M)
  - CON: Regression constant term
  - RESID: Residual sequence (N)
- Regression model equations:
  - Original model: $Y(t) = C_0 + \sum_{j=1}^{M} C_j X_j(t) + \varepsilon(t)$
  - Standardized model: $\frac{Y-\bar{Y}}{\sigma_Y} = \sum_{j=1}^{M} C_{norm,j} \frac{X_j-\bar{X_j}}{\sigma_{X_j}}$
- Least squares solution:
  - Matrix form: $\vec{Y} = \mathbf{X}\vec{C} + \vec{\varepsilon}$
  - Normal equations: $\mathbf{A}\vec{C} = \vec{S}$ where $\mathbf{A} = \mathbf{X}^T\mathbf{X}$, $\vec{S} = \mathbf{X}^T\vec{Y}$
  - Solution vector: $\vec{C} = \mathbf{A}^{-1}\vec{S} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\vec{Y}$
- Statistical calculations:
  - Constant term: $C_0 = \bar{Y} - \sum_{j=1}^{M} C_j \bar{X_j}$
  - Standardized coefficients: $C_{norm,j} = C_j \frac{\sigma_{X_j}}{\sigma_Y}$
  - Residuals: $\text{RESID}(i) = Y(i) - \hat{Y}(i)$
  - Residual variance: $\sigma^2 = \frac{1}{N-M} \sum_{i=1}^{N} \text{RESID}(i)^2$
- Covariance matrix calculation:
  - Coefficient covariance: $\text{COV}(\vec{C}) = \sigma^2 (\mathbf{X}^T\mathbf{X})^{-1}$
  - Standard error: $\text{SE}(C_j) = \sqrt{\text{COV}(j,j)}$
- Numerical method characteristics:
  - Gaussian elimination for matrix inversion
  - Pivot selection to avoid numerical singularity
  - Automatic missing value handling and propagation
- Missing value handling strategy:
  - Strict row-wise processing: exclude entire row if any variable is missing
  - Valid sample counting: use only complete data for regression
  - Missing value flag propagation to residual sequence
- Statistical diagnostic functions:
  - Coefficient of determination: $R^2 = 1 - \frac{\sum \text{RESID}^2}{\sum (Y-\bar{Y})^2}$
  - F-statistic and significance testing support
  - Residual analysis and model validation
- Standardized regression significance:
  - Eliminates dimensional effects, facilitates comparison of relative importance of variables
  - Absolute value of standardized coefficients reflects variable contribution
  - Suitable for multi-variable importance ranking analysis
- Features:
  - Complete implementation of multivariate linear regression
  - Provides both original and standardized coefficients
  - Strict numerical stability guarantees
  - Comprehensive statistical diagnostic output
  - Flexible missing value handling mechanism
  - Suitable for climate index modeling, forecast equation development, and multivariate correlation analysis

## 190. ompgen.F90
**Function**: OpenMP parallel code automatic generator (Fortran90)

**Detailed Description**:
- Functional features:
  - Automatic OpenMP directive insertion
  - Fortran90 free-format source code processing
  - Parallel code pattern generation
  - Template-based code generation support
  - Optimized loop parallelization

## 191. triple2grid.f [📝 src](fortran/triple2grid.f)
**Function**: Triplet scattered data to regular grid interpolation algorithm

**Detailed Description**:
- Contains three subroutines:
  - TRIPLE2GRID1: Main interpolation interface program
  - TRIP2GRD2: Fast nearest neighbor interpolation algorithm
  - TRIP2GRD3: Full search nearest neighbor interpolation algorithm
- Input parameters:
  - KZ: Total number of input data points
  - XI,YI,ZI: Input triplet data (x-coordinate, y-coordinate, z-value)
  - ZMSG: Missing value code
  - MX,NY: Output grid dimensions
  - GX,GY: Output grid coordinate arrays
  - DOMAIN: Domain extension coefficient
  - LOOP: Algorithm selection flag (0=TRIP2GRD2, 1=TRIP2GRD3)
  - METHOD: Distance calculation method (0=Euclidean distance, 1=spherical distance)
  - DISTMX: Maximum search distance
- Output parameters:
  - GRID: Interpolated regular grid (MX×NY)
  - IER: Error code
- Interpolation algorithm principles:
  - Euclidean distance: $d = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}$
  - Spherical distance: $d = R_e \times \arccos(\sin\phi_i \sin\phi_j + \cos\phi_i \cos\phi_j \cos(\lambda_i - \lambda_j))$
  - Nearest neighbor assignment: Assign each observation point to the closest grid point
- Algorithm optimization features:
  - Uniform grid detection: Automatically detect if X/Y directions are uniformly spaced to optimize search
  - Missing value preprocessing: Automatically exclude missing value points
  - External point handling: Handle data points outside grid boundaries through domain extension
  - Exact match priority: Points with exact coordinate matches are processed first
- Distance calculation options:
  - METHOD=0: Cartesian coordinate system, suitable for small-area projected data
  - METHOD=1: Spherical coordinate system, suitable for latitude-longitude data, using Earth radius $R_e = 6371.22$ km
- Boundary handling strategies:
  - In-domain data: Direct interpolation on target grid
  - Out-of-domain data: Create extended grid, then extract target region
  - Extension factor: DOMAIN parameter controls boundary extension range
- Performance comparison:
  - TRIP2GRD2: Fast algorithm, suitable for most cases
  - TRIP2GRD3: Full search algorithm, computationally intensive but better coverage
- Functional features:
  - Supports irregularly distributed scattered data
  - Adaptive grid spacing detection
  - Double precision floating point calculation accuracy
  - Flexible distance calculation modes
  - Intelligent handling of out-of-boundary points
  - Suitable for meteorological station observation data gridding, ocean data interpolation, geographic information systems

## 192. grid2triple.f [📝 src](fortran/grid2triple.f)
**Function**: Regular grid data conversion to triplet (three-column) format

**Detailed Description**:
- Contains one subroutine:
  - GRID2TRIPLE: Converts 2D grid data to triplet list
- Input parameters:
  - X(MX): X-direction coordinate array
  - Y(NY): Y-direction coordinate array
  - Z(MX,NY): 2D data grid
  - MX,NY: Grid dimensions
  - LDMAX: Maximum capacity of output array (typically MX×NY)
  - ZMSG: Missing value code
- Output parameters:
  - D(LDMAX,3): Triplet data array, each row contains [x,y,z] values
  - LD: Actual number of valid data points (excluding missing values)
  - IER: Error code (-10 indicates all values are missing)
- Conversion process:
  - Traverse 2D grid in row-major order: Z(1,1), Z(2,1), ..., Z(MX,1), Z(1,2), ...
  - Skip missing value points: only save data points where Z(m,n) $\neq$ ZMSG
  - Output format: D(i,1)=X(m), D(i,2)=Y(n), D(i,3)=Z(m,n)
- Data structure transformation:
  - Input: Structured grid Z[MX×NY] + coordinate vectors X[MX], Y[NY]
  - Output: Unstructured point cloud [(x₁,y₁,z₁), (x₂,y₂,z₂), ..., (xₙ,yₙ,zₙ)]
- Application scenarios:
  - Grid data preprocessing: Prepare input data for scattered point interpolation algorithms
  - Data format conversion: From regular grid to point cloud format
  - Statistical analysis: Facilitate point-by-point statistics on grid data
  - Visualization: Convert to format suitable for scatter plots
- Functional features:
  - Automatically filters missing values, improving data quality
  - Maintains coordinate-value correspondence
  - Memory-efficient single-pass scanning algorithm
  - Supports arbitrary size regular grids
  - Suitable for climate data analysis, geographic information processing, scientific visualization

## 193. linint2.f [📝 src](fortran/linint2.f)
**Function**: One-dimensional and two-dimensional linear interpolation algorithm collection

**Detailed Description**:
- Contains multiple subroutines:
  - DLININT1: Main program for 1D piecewise linear interpolation
  - DLININT2: Main program for 2D bilinear interpolation
  - DLININT2PTS: 2D point-to-point interpolation program
  - DLIN2INT1: Core algorithm for 1D interpolation
  - DLINT2XY: Core algorithm for 2D interpolation
  - DLINCYC: Cyclic boundary handling
  - DMONOID1/DMONOID2: Monotonicity verification
  - ESTFOW: Missing value estimation algorithm
- One-dimensional interpolation (DLININT1):
  - Input parameters: XI(NXI), FI(NXI), XO(NXO), ICYCX, IOPT
  - Output parameters: FO(NXO), IER
  - Supports cyclic boundaries: ICYCX=0(non-cyclic), $\neq$0(cyclic)
  - Interpolation options: IOPT=0(preserve missing regions), IOPT=1(fill as much as possible)
- Two-dimensional interpolation (DLININT2):
  - Input parameters: XI(NXI), YI(NYI), FI(NXI,NYI), XO(NXO), YO(NYO)
  - Output parameters: FO(NXO,NYO), IER
  - Algorithm flow: First interpolate in X direction → Then interpolate in Y direction
  - Bilinear interpolation formula: $f(x,y) = f(x_1,y_1)(1-t_x)(1-t_y) + f(x_2,y_1)t_x(1-t_y) + f(x_1,y_2)(1-t_x)t_y + f(x_2,y_2)t_x t_y$
  - Where: $t_x = \frac{x-x_1}{x_2-x_1}$, $t_y = \frac{y-y_1}{y_2-y_1}$
- Point-to-point interpolation (DLININT2PTS):
  - Input coordinate pairs: XO(NXYO), YO(NXYO)
  - Output: FO(NXYO)
  - Suitable for irregularly distributed target points
- Cyclic boundary handling characteristics:
  - Automatic boundary extension: Add virtual points at both ends of array
  - Maintain periodicity: Left boundary = right boundary value, right boundary = left boundary value
  - Spacing estimation: Estimate boundary spacing based on existing grid spacing
- Missing value handling strategy:
  - Exact match: Direct assignment when coordinates are identical
  - Complete interpolation: Execute standard bilinear interpolation when all four corner points are valid
  - Partial interpolation: Use distance-weighted averaging when some corner points are missing
  - Distance weights: $w_{i,j} = \frac{1}{\sqrt{(x_i-x_0)^2+(y_j-y_0)^2}}$
- Monotonicity verification:
  - DMONOID1: Verify monotonicity of single sequence (increasing/decreasing)
  - DMONOID2: Verify monotonicity consistency of two sequences
  - Supports strictly monotonic increasing and strictly monotonic decreasing
- Functional features:
  - Efficient piecewise linear interpolation algorithm
  - Complete cyclic boundary support
  - Intelligent missing value handling
  - Strict input validation and error checking
  - Double precision floating point calculation accuracy
  - Suitable for meteorological data interpolation, oceanographic gridding, earth science numerical simulation

## 194. linmsg_dp.f [📝 src](fortran/linmsg_dp.f)
**Function**: Time series missing value linear interpolation filling algorithm

**Detailed Description**:
- Contains one subroutine:
  - DLINMSG: Main program for 1D time series missing value linear interpolation
- Input parameters:
  - X(NPTS): Input time series array (may contain missing values)
  - NPTS: Series length
  - XMSG: Missing value code
  - MFLAG: Endpoint handling flag
    - MFLAG < 0: Set beginning and ending missing values to nearest non-missing value
    - MFLAG $\geq$ 0: Keep beginning and ending missing values as missing
  - MPTCRT: Maximum consecutive missing value threshold
    - When consecutive missing values exceed this threshold, no interpolation is performed
    - Usually set to NPTS to interpolate as much as possible
- Output results:
  - X(NPTS): Interpolated time series (in-place modification)
  - Return status code tracked through internal variables
- Algorithm flow:
  1. **Complete missing check**: If entire series is missing values, return directly
  2. **Missing segment identification**: Scan series to mark start (NSTRT) and end (NEND) positions of consecutive missing value segments
  3. **Threshold check**: If missing segment length > MPTCRT, skip that segment without interpolation
  4. **Beginning interpolation**: If series starts with missing values (NSTRT=1)
     - MFLAG < 0: Fill with first valid value
     - MFLAG $\geq$ 0: Keep missing state
  5. **Middle segment interpolation**: Use linear interpolation formula
     - Slope: $\text{slope} = \frac{X(N) - X(\text{NBASE})}{N - \text{NBASE}}$
     - Interpolation: $X(NN) = X(\text{NBASE}) + \text{slope} \times (NN - \text{NBASE})$
  6. **End interpolation**: If series ends with missing values (NEND=NPTS)
     - MFLAG < 0: Fill with last valid value
     - MFLAG $\geq$ 0: Keep missing state
- Interpolation characteristics:
  - Maintains temporal continuity of data
  - Linear trend preservation: interpolated points lie strictly on the line connecting valid points
  - Flexible boundary handling: can choose to fill or maintain missing state
  - Threshold control: avoid unreasonable interpolation across excessively long missing segments
- Application scenarios:
  - Meteorological observation data preprocessing: fill data gaps during instrument failures
  - Time series analysis: prepare continuous data for statistical analysis
  - Data quality control: reasonable range missing value repair
  - Model input preparation: provide complete time series for numerical models
- Functional features:
  - In-place operation, memory efficient
  - Supports time series of arbitrary length
  - Flexible endpoint handling strategy
  - Configurable interpolation threshold control
  - Double precision floating point calculation accuracy
  - Suitable for preprocessing and quality control of various time series data

## 195. moc_loops.f [📝 src](fortran/moc_loops.f)
**Function**: Meridional Overturning Circulation (MOC) regional integral calculation

**Detailed Description**:
- Contains one subroutine:
  - MOCLOOPS: Main program for latitude band and depth layer integration of ocean meridional overturning circulation
- Input parameters:
  - NYAUX: Number of auxiliary latitude grid points
  - MLON, NLAT, KDEP: Longitude, latitude, and depth dimensions
  - NRX: Number of regions (1=global, 2=Atlantic)
  - TLAT(MLON, NLAT): Actual latitude values at each grid point
  - LAT_AUX_GRID(NYAUX): Auxiliary latitude grid boundaries
  - RMLAK(MLON, NLAT, 2): Regional mask array
    - RMLAK(:,:,1)=1: Global ocean points
    - RMLAK(:,:,2)=1: Atlantic ocean points
  - WORK1, WORK2, WORK3(MLON, NLAT, KDEP): Three input 3D ocean data fields
  - WMSG: Missing value code
- Output parameters:
  - TMP1, TMP2, TMP3(NYAUX, KDEP, 2): Latitude-depth-region integral results
    - First dimension: Latitude band index
    - Second dimension: Depth layer index
    - Third dimension: Region index (1=global, 2=Atlantic)
- Algorithm flow:
  1. **Initialization**: Set all output arrays TMP1, TMP2, TMP3 to 0
  2. **Global integration** (NR=1):
     - Traverse latitude bands: LAT_AUX_GRID(NY-1) $\leq$ TLAT < LAT_AUX_GRID(NY)
     - Process only global ocean points: RMLAK(ML, NL, 1)=1
     - Skip missing values: WORK1(ML, NL, KD) $\neq$ WMSG
     - Accumulate sum: TMP1(NY, KD, 1) += WORK1(ML, NL, KD)
  3. **Atlantic integration** (NR=2):
     - Same latitude band grouping logic
     - Process only Atlantic ocean points: RMLAK(ML, NL, 2)=1
     - Accumulate results to TMP1(NY, KD, 2), TMP2(NY, KD, 2), TMP3(NY, KD, 2)
- MOC calculation principle:
  - Meridional overturning circulation reflects north-south mass transport in the ocean
  - Integration by latitude band removes zonal (east-west) variation
  - Separation by depth layer distinguishes circulation structure at different depths
  - Regional separation enables comparison of global and basin-scale circulation differences
- Data processing features:
  - Automatic latitude band grouping: grid points are classified based on actual latitude and auxiliary grid
  - Regional selective integration: mask array enables land-sea and basin separation
  - Safe missing value handling: invalid data points are automatically skipped
  - Parallel calculation of three fields: can process velocity, temperature, salinity, etc. simultaneously
- Application scenarios:
  - Climate model post-processing: calculate MOC strength from ocean model output
  - Ocean circulation analysis: study water exchange at different latitudes and depths
  - Climate change monitoring: track long-term circulation trends
  - Model comparison and validation: evaluate MOC performance in different ocean models
- Functional features:
  - Efficient triple-nested loop algorithm
  - Flexible regional mask support
  - Simultaneous multi-field processing
  - Automatic latitude band grouping mechanism
  - Double precision floating point accuracy
  - Suitable for MOC diagnostics in global ocean models

## 196. rcm2points.f [📝 src](fortran/rcm2points.f)
**Function**: Interpolation algorithm from regional climate model grid to scattered points

**Detailed Description**:
- Contains one subroutine:
  - DRCM2POINTS: Interpolates regular/irregular 2D grid data to specified scattered point locations
- Input parameters:
  - NGRD: Number of data layers (can process multiple variables simultaneously)
  - NXI, NYI: Input grid dimensions
  - XI(NXI, NYI), YI(NXI, NYI): Input grid longitude and latitude coordinates (can be irregular grid)
  - FI(NXI, NYI, NGRD): Multi-layer input data fields
  - NXYO: Number of output scattered points
  - XO(NXYO), YO(NXYO): Target scattered point longitude and latitude coordinates
  - XMSG: Missing value code
  - OPT: Interpolation method option (0/1 = inverse distance weighting, 2 = bilinear interpolation)
  - NCRIT: Valid point threshold (usually ≥3)
  - KVAL: Grid step parameter (usually = 1)
- Output parameters:
  - FO(NXYO, NGRD): Interpolated multi-layer scattered point data
  - IER: Error code
- Three-stage interpolation algorithm:
  1. **Exact match stage**:
     - Find points with exactly matching coordinates: XO = XI and YO = YI
     - Direct assignment: FO(NXY, NG) = FI(IX, IY, NG)
  2. **Local interpolation stage**:
     - Search for target point within a 2×2 grid cell
     - Condition: XI(IX, IY) ≤ XO ≤ XI(IX+K, IY) and YI(IX, IY) ≤ YO ≤ YI(IX, IY+K)
     - Bilinear interpolation (OPT=2): $w_{i,j} = (1-t_x)(1-t_y), t_x(1-t_y), (1-t_x)t_y, t_x t_y$
     - Inverse distance weighting (OPT≠2): $w_{i,j} = \frac{1}{d_{i,j}^2}$, where $d_{i,j}$ is spherical distance
  3. **Global search stage**:
     - For points still not interpolated, perform radius search
     - Search radius: DKM = 5° latitude ≈ 555 km
     - Condition: |YI - YO| ≤ 5° and spherical distance ≤ 555 km
     - Inverse distance weighting: $w = \frac{1}{d^2}$, $FO = \frac{\sum w_i \cdot FI_i}{\sum w_i}$
- Distance calculation:
  - Uses DGCDIST function for spherical great-circle distance
  - Formula: $d = R_e \times \arccos(\sin\phi_1 \sin\phi_2 + \cos\phi_1 \cos\phi_2 \cos(\lambda_1-\lambda_2))$
  - Earth radius: RE = 6371 km
- Quality control features:
  - Monotonicity check: Verifies monotonicity of input grid longitude and latitude
  - Valid point check: Requires at least NCRIT valid neighboring points for interpolation
  - Missing value handling: Automatically skips missing values and adjusts weights
  - Multi-layer synchronization: All data layers use the same interpolation weights
- Application scenarios:
  - Regional climate model post-processing: Interpolate model grid data to observation stations
  - Model-observation comparison: Extract model grid values at observation point locations
  - Data assimilation: Prepare observation operators for assimilation systems
  - Sensitivity analysis: Extract model output at specific locations
- Functional features:
  - Supports irregular grid input (curvilinear grid)
  - Multiple interpolation methods available (bilinear/inverse distance)
  - Three-stage progressive interpolation strategy
  - Simultaneous processing of multiple data fields
  - Strict quality control and error checking
  - Suitable for scattered point interpolation in regional climate, ocean, and atmospheric models

## 197. rcm2rgrid.f [📝 src](fortran/rcm2rgrid.f)
**Function**: Bidirectional interpolation algorithm between regional climate model (RCM) grids and regular grids

**Detailed Description**:
- Contains three subroutines:
  - DRCM2RGRID: Interpolation from curvilinear grid to regular grid
  - DRGRID2RCM: Interpolation from regular grid to curvilinear grid  
  - DGCDIST: Spherical great-circle distance calculation function
- DRCM2RGRID interpolation features:
  - Input: Curvilinear grid XI(NXI,NYI), YI(NXI,NYI), FI(NXI,NYI,NGRD)
  - Output: Regular grid XO(NXO), YO(NYO), FO(NXO,NYO,NGRD)
  - Three-stage interpolation strategy:
    1. **Approximate exact match**: Tolerance EPS=1×10⁻⁴, finds nearly overlapping points
    2. **Inverse distance weighting**: Uses spherical inverse distance weights within 2×2 curvilinear grid cells
    3. **Linear interpolation fill**: Performs linear interpolation along longitude for remaining missing points
- DRGRID2RCM interpolation features:
  - Input: Regular grid XI(NXI), YI(NYI), FI(NXI,NYI,NGRD)
  - Output: Curvilinear grid XO(NXO,NYO), YO(NXO,NYO), FO(NXO,NYO,NGRD)
  - Two-stage interpolation strategy:
    1. **Approximate exact match**: Tolerance EPS=1×10⁻³
    2. **Bilinear/inverse distance weighting**: 
       - All four corner points valid: standard bilinear interpolation FBLI()
       - Partial missing: inverse distance weighting interpolation
- Bilinear interpolation formula:
  - Inline function: FLI(Z1,Z2,SLOPE) = Z1 + SLOPE×(Z2-Z1)
  - 2D extension: FBLI(Z1,Z2,Z3,Z4,SLPX,SLPY) = interpolate in X, then in Y
  - Weight calculation: SLPX = (XO-XI₁)/(XI₂-XI₁), SLPY = (YO-YI₁)/(YI₂-YI₁)
- DGCDIST spherical distance algorithm:
  - Input: Two points' lat/lon (RLAT1,RLON1), (RLAT2,RLON2)
  - Output unit selection: IU=1(radian), IU=2(degree), IU=3(meter), IU=4(kilometer)
  - Calculation formula: Haversine variant, handles longitude crossing 180°
  - Special optimization: Identical point detection returns zero distance
  - Earth radius: 6371.22 km
- Interpolation algorithm optimizations:
  - Monotonicity pre-check: Ensures monotonicity of input/output grid coordinates
  - Exact match priority: Avoids unnecessary numerical interpolation errors
  - Threshold control: NCRIT parameter controls minimum number of valid neighbors
  - Missing value safety: Complete missing value checks and weight normalization
- Curvilinear grid handling features:
  - Supports arbitrary curvilinear coordinate systems (e.g., Lambert projection, spherical grids)
  - Automatically detects grid cell containment
  - Handles grid distortion and irregular spacing
  - Three-stage progressive filling strategy ensures high coverage
- Application scenarios:
  - RCM post-processing: Remap RCM output to standard lat-lon grids
  - Model nesting: Interpolate global model data to regional model grids
  - Observation assimilation: Regrid satellite observation data
  - Data standardization: Unify grids between different models
- Functional features:
  - Bidirectional interpolation capability (curvilinear ↔ regular)
  - Adaptive selection of multiple interpolation methods
  - High-precision spherical geometry calculations
  - Batch multi-layer data processing
  - Robust quality control mechanism
  - Suitable for numerical weather prediction, climate modeling, oceanographic research


# Summary

This document covers 197 Fortran source files in the NCL project, including:
- 173 Fortran 77 format files (.f)
- 21 Fortran 90/95 format files (.f90)
- 1 Fortran 90 module file (.F90)
- 1 Python interface definition file (.pyf)

These files span scientific computing algorithms in meteorology, oceanography, climatology, numerical analysis, statistics, and other fields.

The files cover the following major functional areas:

## Core Functional Categories:

### 1. Meteorology & Atmospheric Science (30%)
- Temperature, humidity, and pressure conversion
- Atmospheric physical process calculations
- WRF model-specific functions
- Vertical coordinate transformations
- Atmospheric dynamics computations

### 2. Numerical Analysis & Computation (25%)
- FFT transforms and frequency domain analysis
- Interpolation and extrapolation algorithms
- Numerical integration and differentiation
- Matrix operations and linear algebra
- Solvers and optimization algorithms

### 3. Statistical Analysis (20%)
- Probability distribution functions
- Regression analysis and correlation
- Time series analysis
- Hypothesis testing
- Clustering and classification

### 4. Data Processing & Gridding (15%)
- Grid conversion and remapping
- Missing value handling
- Data averaging and aggregation
- Quality control
- File I/O operations

### 5. Spherical Geometry & Earth Science (10%)
- Spherical harmonic analysis
- Geographic coordinate transformations
- Great circle calculations
- Geophysical parameters
#   M e t e o r o l o g y - F o r t r a n - F u n c t i o n s  
 