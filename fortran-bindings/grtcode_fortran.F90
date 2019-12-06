!> @file
!GRTCODE is a GPU-able Radiative Transfer Code
!Copyright (C) 2016  Garrett Wright
!Modified in 2019 by Raymond Menzel
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; version 2.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

!> @brief Fortran bindings for utilities.
module grtcode
use, intrinsic :: iso_c_binding, only: c_char, c_double, c_float, c_int, c_int64_t, &
                                       c_null_char, c_null_ptr, c_ptr
implicit none
private


#ifdef SINGLE_PRECISION
integer, parameter :: fp = c_float
#else
integer, parameter :: fp = c_double
#endif
integer, parameter, public :: grtcode_success = 0
integer, parameter :: grid_struct = 0
integer, parameter :: optics_struct = 1
integer, parameter :: gas_optics_struct = 2
integer, parameter :: solar_flux_struct = 3
integer(kind=c_int), parameter, public :: H2O = 1
integer(kind=c_int), parameter, public :: CO2 = 2
integer(kind=c_int), parameter, public :: O3 = 3
integer(kind=c_int), parameter, public :: N2O = 4
integer(kind=c_int), parameter, public :: CO = 5
integer(kind=c_int), parameter, public :: CH4 = 6
integer(kind=c_int), parameter, public :: O2 = 7
integer(kind=c_int), parameter, public :: NO = 8
integer(kind=c_int), parameter, public :: SO2 = 9
integer(kind=c_int), parameter, public :: NO2 = 10
integer(kind=c_int), parameter, public :: NH3 = 11
integer(kind=c_int), parameter, public :: HNO3 = 12
integer(kind=c_int), parameter, public :: OH = 13
integer(kind=c_int), parameter, public :: HF = 14
integer(kind=c_int), parameter, public :: HCl = 15
integer(kind=c_int), parameter, public :: HBr = 16
integer(kind=c_int), parameter, public :: HI = 17
integer(kind=c_int), parameter, public :: ClO = 18
integer(kind=c_int), parameter, public :: OCS = 19
integer(kind=c_int), parameter, public :: H2CO = 20
integer(kind=c_int), parameter, public :: HOCl = 21
integer(kind=c_int), parameter, public :: N2 = 22
integer(kind=c_int), parameter, public :: HCN = 23
integer(kind=c_int), parameter, public :: CH3Cl = 24
integer(kind=c_int), parameter, public :: H2O2 = 25
integer(kind=c_int), parameter, public :: C2H2 = 26
integer(kind=c_int), parameter, public :: C2H6 = 27
integer(kind=c_int), parameter, public :: PH3 = 28
integer(kind=c_int), parameter, public :: COF2 = 29
integer(kind=c_int), parameter, public :: SF6_MOL = 30
integer(kind=c_int), parameter, public :: H2S = 31
integer(kind=c_int), parameter, public :: HCOOH = 32
integer(kind=c_int), parameter, public :: HO2 = 33
integer(kind=c_int), parameter, public :: O = 34
integer(kind=c_int), parameter, public :: ClONO2 = 35
integer(kind=c_int), parameter, public :: NOp = 36
integer(kind=c_int), parameter, public :: HOBr = 37
integer(kind=c_int), parameter, public :: C2H4 = 38
integer(kind=c_int), parameter, public :: CH3OH = 39
integer(kind=c_int), parameter, public :: CH3Br = 40
integer(kind=c_int), parameter, public :: CH3CN = 41
integer(kind=c_int), parameter, public :: CF4_MOL = 42
integer(kind=c_int), parameter, public :: C4H2 = 43
integer(kind=c_int), parameter, public :: HC3N = 44
integer(kind=c_int), parameter, public :: H2 = 45
integer(kind=c_int), parameter, public :: CS = 46
integer(kind=c_int), parameter, public :: SO3 = 47
integer(kind=c_int), parameter, public :: C2N2 = 48
integer(kind=c_int), parameter, public :: COCl2 = 49
integer(kind=c_int), parameter, public :: SO = 50
integer(kind=c_int), parameter, public :: C3H4 = 51
integer(kind=c_int), parameter, public :: CH3 = 52
integer(kind=c_int), parameter, public :: CS2 = 53
integer(kind=c_int), parameter, public :: MAX_NUM_MOLECULES = 53
integer(kind=c_int), parameter, public :: CFC11 = 0
integer(kind=c_int), parameter, public :: CFC12 = 1
integer(kind=c_int), parameter, public :: CFC113 = 2
integer(kind=c_int), parameter, public :: CFC114 = 3
integer(kind=c_int), parameter, public :: CFC115 = 4
integer(kind=c_int), parameter, public :: HCFC22 = 5
integer(kind=c_int), parameter, public :: HCFC141b = 6
integer(kind=c_int), parameter, public :: HCFC142b = 7
integer(kind=c_int), parameter, public :: HFC23 = 8
integer(kind=c_int), parameter, public :: HFC125 = 9
integer(kind=c_int), parameter, public :: HFC134a = 10
integer(kind=c_int), parameter, public :: HFC143a = 11
integer(kind=c_int), parameter, public :: HFC152a = 12
integer(kind=c_int), parameter, public :: HFC227ea = 13
integer(kind=c_int), parameter, public :: HFC245fa = 14
integer(kind=c_int), parameter, public :: CCl4 = 15
integer(kind=c_int), parameter, public :: C2F6 = 16
integer(kind=c_int), parameter, public :: CF4 = 17
integer(kind=c_int), parameter, public :: CH2Cl2 = 18
integer(kind=c_int), parameter, public :: NF3 = 19
integer(kind=c_int), parameter, public :: SF6 = 20
integer(kind=c_int), parameter, public :: MAX_NUM_CFCS = 21
integer(kind=c_int), parameter, public :: CIA_N2 = 0
integer(kind=c_int), parameter, public :: CIA_O2 = 1
integer(kind=c_int), parameter, public :: MAX_NUM_CIAS = 2


!> @brief Device object.
type, public :: Device_t
  integer(kind=c_int) :: device !< Device identifier.
end type Device_t


interface create_device
  function c_create_device(device, id) &
    result(error_code) &
    bind(c, name="create_device")
    import c_int
    integer(kind=c_int), intent(inout) :: device
    integer(kind=c_int), intent(in), optional :: id
    integer(kind=c_int) :: error_code
  end function c_create_device
  module procedure f_create_device
end interface create_device
public :: create_device


!> @brief One dimensional grid.
type, public :: Grid_t
  type(c_ptr) :: grid !< Pointer to grid object.
end type Grid_t


interface create_spectral_grid
  function c_create_spectral_grid(grid, w0, wn, dw) &
    result(error_code) &
    bind(c, name="create_spectral_grid")
    import c_double, c_int, c_ptr
    type(c_ptr), value :: grid
    real(kind=c_double), intent(in), value :: w0
    real(kind=c_double), intent(in), value :: wn
    real(kind=c_double), intent(in), value :: dw
    integer(kind=c_int) :: error_code
  end function c_create_spectral_grid
  module procedure f_create_spectral_grid
end interface create_spectral_grid
public :: create_spectral_grid
public :: destroy_spectral_grid


!> @brief Gas optical properties.
type, public :: Optics_t
  type(c_ptr) :: optics !< Pointer to optics object.
end type Optics_t


interface create_optics
  function c_create_optics(optics, num_layers, grid, device) &
    result(error_code) &
    bind(c, name="create_optics")
    import c_int, c_ptr
    type(c_ptr), value :: optics
    integer(kind=c_int), intent(in), value :: num_layers
    type(c_ptr), intent(in), value :: grid
    integer(kind=c_int), intent(in) :: device
    integer(kind=c_int) :: error_code
  end function c_create_optics
  module procedure f_create_optics
end interface create_optics
public :: create_optics


interface destroy_optics
  function c_destroy_optics(optics) &
    result(error_code) &
    bind(c, name="destroy_optics")
    import c_int, c_ptr
    type(c_ptr), value :: optics
    integer(kind=c_int) :: error_code
  end function c_destroy_optics
  module procedure f_destroy_optics
end interface destroy_optics
public :: destroy_optics


interface
  function malloc_struct(p, type_) &
    result(error_code) &
    bind(c)
    import c_int, c_ptr
    type(c_ptr) :: p
    integer(kind=c_int), intent(in), value :: type_
    integer(kind=c_int) :: error_code
  end function malloc_struct
end interface
public :: malloc_struct


interface
  function free_struct(p) &
    result(error_code) &
    bind(c)
    import c_int, c_ptr
    type(c_ptr) :: p
    integer(kind=c_int) :: error_code
  end function free_struct
end interface
public :: free_struct


interface
  subroutine grtcode_set_verbosity(level) &
    bind(c)
    import c_int
    integer(kind=c_int), intent(in), value :: level
  end subroutine grtcode_set_verbosity
end interface
public :: grtcode_set_verbosity


interface optical_properties
  function c_optical_properties(optics, tau, omega, g) &
    result(error_code) &
    bind(c, name="optical_properties")
    import c_int, c_ptr, fp
    type(c_ptr), intent(in), value :: optics
    real(kind=fp), dimension(*), intent(inout), optional :: tau
    real(kind=fp), dimension(*), intent(inout), optional :: omega
    real(kind=fp), dimension(*), intent(inout), optional :: g
    integer(kind=c_int) :: error_code
  end function c_optical_properties
  module procedure f_optical_properties
end interface optical_properties
public :: optical_properties


interface spectral_grid_properties
  !> @brief Get the spectral grid properties.
  !! @return RS_SUCCESS or an error code.
  function c_spectral_grid_properties(grid, w0, n, dw) &
    result(return_code) &
    bind(c, name="spectral_grid_properties")
    import c_double, c_int, c_int64_t, c_ptr
    type(c_ptr), intent(in), value :: grid !< Spectral grid.
    real(kind=c_double), intent(out), optional :: w0 !< Grid lower bound.
    integer(kind=c_int64_t), intent(out), optional :: n !< Grid size.
    real(kind=c_double), intent(out), optional :: dw !< Grid spacing.
    integer(kind=c_int) :: return_code
  end function c_spectral_grid_properties
  module procedure f_spectral_grid_properties
end interface spectral_grid_properties
public :: spectral_grid_properties


interface add_optics
  !> @brief Add optical properties together.
  !! @return RS_SUCCESS or an error code.
  function c_add_optics(optics, num_optics, res) &
    result(return_code) &
    bind(c, name="add_optics")
    import c_int, c_ptr
    type(c_ptr), dimension(*) :: optics
    integer(kind=c_int), intent(in), value :: num_optics
    type(c_ptr), value :: res
    integer(kind=c_int) :: return_code
  end function c_add_optics
  module procedure f_add_optics
end interface add_optics
public :: add_optics


type, public :: SolarFlux_t
  type(c_ptr) :: solar_flux !< Pointer to solar flux object.
end type SolarFlux_t


interface create_solar_flux
  !> @brief Read in data for the solar flux.
  !! @return RS_SUCCESS or an error code.
  function c_create_solar_flux(solar_flux, grid, path) &
    result(return_code) &
    bind(c, name="create_solar_flux")
    import c_char, c_int, c_ptr
    type(c_ptr), value :: solar_flux !< Solar flux object.
    type(c_ptr), value :: grid !< Spectral grid.
    character(kind=c_char, len=1), dimension(*), intent(in) :: path !< Solar flux csv file.
    integer(kind=c_int) :: return_code
  end function c_create_solar_flux
  module procedure f_create_solar_flux
end interface create_solar_flux
public :: create_solar_flux


interface destroy_solar_flux
  !> @brief Free memory for the solar flux.
  !! @return RS_SUCCESS or an error code.
  function c_destroy_solar_flux(solar_flux) &
    result(return_code) &
    bind(c, name="destroy_solar_flux")
    import c_int, c_ptr
    type(c_ptr), value :: solar_flux !< Solar flux object.
    integer(kind=c_int) :: return_code
  end function c_destroy_solar_flux
  module procedure f_destroy_solar_flux
end interface destroy_solar_flux
public :: destroy_solar_flux


interface solar_flux_properties
  !> @brief Get the solar flux properties.
  !! @return RS_SUCCESS or an error code.
  function c_solar_flux_properties(solar, flux) &
    result(return_code) &
    bind(c, name="solar_flux_properties")
    import c_int, c_ptr, fp
    type(c_ptr), intent(in), value :: solar !< Solar flux.
    real(kind=fp), dimension(*), intent(inout) :: flux !< Flux.
    integer(kind=c_int) :: return_code
  end function c_solar_flux_properties
  module procedure f_solar_flux_properties
end interface solar_flux_properties
public :: solar_flux_properties


interface get_num_gpus
  !> @brief Determine the number of CUDA-enabled GPUs on the system.
  !! @return RS_SUCCESS or an error code.
  function c_get_num_gpus(num_devices, verbose) &
    result(return_code) &
    bind(c, name="get_num_gpus")
    import c_int
    integer(kind=c_int), intent(inout) :: num_devices !< Number of CUDA-enabled devices found.
    integer(kind=c_int), intent(in), value :: verbose !< Verbosity flag.
    integer(kind=c_int) :: return_code
  end function c_get_num_gpus
end interface get_num_gpus
public :: get_num_gpus


type, public :: MolecularLines_t
  type(c_ptr) :: ml !< Pointer to molecular lines object.
end type MolecularLines_t


interface create_gas_optics
  !> @brief Reserve memory for molecular lines.
  !! @return RS_SUCCESS or an error code.
  function c_create_gas_optics(ml, num_levels, grid, device, hitran_path, h2o_ctm_dir, &
                               o3_ctm_file, wcutoff, optical_depth_method) &
    result(return_code) &
    bind(c, name="create_gas_optics")
    import c_char, c_double, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: num_levels !< Number of atmospheric levels.
    type(c_ptr), intent(in), value :: grid !< Spectral grid.
    integer(kind=c_int), intent(in) :: device !< Device.
    character(kind=c_char, len=1), dimension(*), intent(in) :: hitran_path !< Path to HITRAN database file.
    character(kind=c_char, len=1), dimension(*), intent(in) :: h2o_ctm_dir !< Path to water vapor continuum directory.
    character(kind=c_char, len=1), dimension(*), intent(in) :: o3_ctm_file !< Path to ozone continuum file.
    real(kind=c_double), intent(in), optional :: wcutoff !< Cutoff from line center [1/cm].
    integer(kind=c_int), intent(in), optional :: optical_depth_method !< Method used to calculate the optical depths.
    integer(kind=c_int) :: return_code
  end function c_create_gas_optics
  module procedure f_create_gas_optics
end interface create_gas_optics
public :: create_gas_optics


interface destroy_gas_optics
  !> @brief Free memory for the molecular lines.
  !! @return RS_SUCCESS or an error code.
  function c_destroy_gas_optics(ml) &
    result(return_code) &
    bind(c, name="destroy_gas_optics")
    import c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int) :: return_code
  end function c_destroy_gas_optics
  module procedure f_destroy_gas_optics
end interface destroy_gas_optics
public :: destroy_gas_optics


interface add_molecule
  !> @brief Add a molecule.
  !! @return RS_SUCCESS or an error code.
  function grt_add_molecule(ml, molecule_id, min_line_center, max_line_center) &
    result(return_code) &
    bind(c, name="add_molecule")
    import c_double, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: molecule_id !< Molecule id.
    real(kind=c_double), intent(in), optional :: min_line_center !< Lower bound [1/cm] for spectral line centers.
    real(kind=c_double), intent(in), optional :: max_line_center !< Upper bound [1/cm] for spectral line centers.
    integer(kind=c_int) :: return_code
  end function grt_add_molecule
  module procedure f_add_molecule
end interface add_molecule
public :: add_molecule


interface set_molecule_ppmv
  !> @brief Update a molecule's ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_molecule_ppmv(ml, molecule_id, ppmv) &
    result(return_code) &
    bind(c, name="set_molecule_ppmv")
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: molecule_id  !< Molecule id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_molecule_ppmv
  module procedure f_set_molecule_ppmv
end interface set_molecule_ppmv
public :: set_molecule_ppmv


interface add_cfc
  !> @brief Add a CFC.
  !! @return RS_SUCCESS or an error code.
  function grt_add_cfc(ml, cfc_id, filepath) &
    result(return_code) &
    bind(c, name="add_cfc")
    import c_char, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: cfc_id !< CFC id.
    character(kind=c_char, len=1), dimension(*), intent(in) :: filepath !< Path to CFC cross section csv file.
    integer(kind=c_int) :: return_code
  end function grt_add_cfc
  module procedure f_add_cfc
end interface add_cfc
public :: add_cfc


interface set_cfc_ppmv
  !> @brief Update a CFC's ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_cfc_ppmv(ml, cfc_id, ppmv) &
    result(return_code) &
    bind(c, name="set_cfc_ppmv")
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: cfc_id !< CFC id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_cfc_ppmv
  module procedure f_set_cfc_ppmv
end interface set_cfc_ppmv
public :: set_cfc_ppmv


interface add_cia
  !> @brief Activate collision-induced absorption between two species.
  !! @return RS_SUCCESS or an error code.
  function grt_add_cia(ml, species1, species2, filepath) &
    result(return_code) &
    bind(c, name="add_cia")
    import c_char, c_int, c_ptr
    type(c_ptr), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(in), value :: species1 !< Id of species.
    integer(kind=c_int), intent(in), value :: species2 !< Id of species.
    character(kind=c_char, len=1), dimension(*), intent(in) :: filepath !< Path to cross section csv file.
    integer(kind=c_int) :: return_code
  end function grt_add_cia
  module procedure f_add_cia
end interface add_cia
public :: add_cia


interface set_cia_ppmv
  !> @brief Update a CIA species' ppmv.
  !! @return RS_SUCCESS or an error code.
  function grt_set_cia_ppmv(ml, cia_id, ppmv) &
    result(return_code) &
    bind(c, name="set_cia_ppmv")
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecularlines object.
    integer(kind=c_int), intent(in), value :: cia_id !< CIA species id.
    real(kind=fp), dimension(*), intent(in) :: ppmv !< Abundance [ppmv] (level).
    integer(kind=c_int) :: return_code
  end function grt_set_cia_ppmv
  module procedure f_set_cia_ppmv
end interface set_cia_ppmv
public :: set_cia_ppmv


interface calculate_optics
  !> @brief Calcluate the total optical depth in each layer at each spectral grid point.
  !! @return RS_SUCCESS or an error code.
  function grt_calculate_optical_depth(ml, pressure, temperature, optics) &
    result(return_code) &
    bind(c, name="calculate_optical_depth")
    import c_int, c_ptr, fp
    type(c_ptr), value :: ml !< Molecular lines object.
    real(kind=fp), dimension(*), intent(in) :: pressure !< Pressure [mb] (level).
    real(kind=fp), dimension(*), intent(in) :: temperature !< Temperature [K] (level).
    type(c_ptr), value :: optics !< Optics object.
    integer(kind=c_int) :: return_code
  end function grt_calculate_optical_depth
  module procedure f_calculate_optics
end interface calculate_optics
public :: calculate_optics


interface num_molecules
  !> @brief Get the number of molecules.
  !! @return RS_SUCCESS or an error code.
  function grt_get_num_molecules(ml, n) &
    result(return_code) &
    bind(c, name="get_num_molecules")
    import c_int, c_ptr
    type(c_ptr), intent(in), value :: ml !< Molecular lines object.
    integer(kind=c_int), intent(out) :: n !< Number of molecules.
    integer(kind=c_int) :: return_code
  end function grt_get_num_molecules
  module procedure f_num_molecules
end interface num_molecules
public :: num_molecules


interface grtcode_errstr
  !> @brief Return a message for an input return code.
  !! @return RS_SUCCESS or an error code.
  function c_grt_errstr(code, buf, buf_size) &
    result(return_code) &
    bind(c, name="grtcode_errstr")
    import c_char, c_int
    integer(kind=c_int), intent(in), value :: code !< Error code.
    character(kind=c_char, len=1), dimension(*) :: buf !< Buffer to hold error message.
    integer(kind=c_int), intent(in), value :: buf_size !< Size of input buffer.
    integer(kind=c_int) :: return_code
  end function c_grt_errstr
  module procedure f_grt_errstr
end interface grtcode_errstr
public :: grtcode_errstr


interface rayleigh_scattering
  !> @brief Calculate the optical properties due to Rayleigh scattering.
  !! @return RS_SUCCESS or an error code.
  function c_rayleigh_scattering(optics, pressure) &
    result(return_code) &
    bind(c, name="rayleigh_scattering")
    import c_int, c_ptr, fp
    type(c_ptr), value :: optics !< Optics object.
    real(kind=fp), dimension(*), intent(in) :: pressure !< Pressure [mb] (level).
    integer(kind=c_int) :: return_code
  end function c_rayleigh_scattering
  module procedure f_rayleigh_scattering
end interface rayleigh_scattering
public :: rayleigh_scattering


contains


subroutine append_null_char(str_in, array_out)
  character(kind=c_char, len=*), intent(in) :: str_in
  character(kind=c_char, len=1), dimension(:), allocatable, intent(inout) :: array_out
  integer :: i
  integer :: s
  if (allocated(array_out)) then
    deallocate(array_out)
  endif
  s = len_trim(str_in)
  allocate(array_out(s+1))
  do i = 1, s
    array_out(i) = str_in(i:i)
  enddo
  array_out(i) = c_null_char
end subroutine append_null_char


function f_create_device(device, id) &
  result(error_code)
  type(Device_t), intent(inout) :: device
  integer(kind=c_int), intent(in), optional :: id
  integer(kind=c_int) :: error_code
  error_code = c_create_device(device%device, id)
end function f_create_device


function f_create_spectral_grid(grid, w0, wn, dw) &
  result(error_code)
  type(Grid_t), intent(inout) :: grid
  real(kind=c_double), intent(in) :: w0
  real(kind=c_double), intent(in) :: wn
  real(kind=c_double), intent(in) :: dw
  integer(kind=c_int) :: error_code
  grid%grid = c_null_ptr
  error_code = malloc_struct(grid%grid, grid_struct)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = c_create_spectral_grid(grid%grid, w0, wn, dw)
end function f_create_spectral_grid


function destroy_spectral_grid(grid) &
  result(error_code)
  type(Grid_t), intent(inout) :: grid
  integer(kind=c_int) :: error_code
  error_code = free_struct(grid%grid)
end function destroy_spectral_grid


function f_create_optics(optics, num_layers, grid, device) &
  result(error_code)
  type(Optics_t), intent(inout) :: optics
  integer(kind=c_int), intent(in) :: num_layers
  type(Grid_t), intent(in) :: grid
  type(Device_t), intent(in) :: device
  integer(kind=c_int) :: error_code
  optics%optics = c_null_ptr
  error_code = malloc_struct(optics%optics, optics_struct)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = c_create_optics(optics%optics, num_layers, grid%grid, device%device)
end function f_create_optics


function f_destroy_optics(optics) &
  result(error_code)
  type(Optics_t), intent(inout) :: optics
  integer(kind=c_int) :: error_code
  error_code = c_destroy_optics(optics%optics)
  if (error_code .ne. grtcode_success) then
    return
  endif
  error_code = free_struct(optics%optics)
end function f_destroy_optics


function f_optical_properties(optics, tau, omega, g) &
  result(error_code)
  type(Optics_t), intent(in) :: optics
  real(kind=fp), dimension(:,:), intent(inout), optional :: tau
  real(kind=fp), dimension(:,:), intent(inout), optional :: omega
  real(kind=fp), dimension(:,:), intent(inout), optional :: g
  integer(kind=c_int) :: error_code
  error_code = c_optical_properties(optics%optics, tau, omega, g)
end function f_optical_properties


function f_spectral_grid_properties(grid, w0, n, dw) &
  result(return_code)
  type(Grid_t), intent(in) :: grid !< Molecular lines object.
  real(kind=c_double), intent(out), optional :: w0 !< Grid lower bound.
  integer(kind=c_int64_t), intent(out), optional :: n !< Grid size.
  real(kind=c_double), intent(out), optional :: dw !< Grid spacing.
  integer(kind=c_int) :: return_code
  return_code = c_spectral_grid_properties(grid%grid, w0, n, dw)
end function f_spectral_grid_properties


function f_add_optics(optics, res) &
  result(return_code)
  type(Optics_t), dimension(:), intent(in) :: optics
  type(Optics_t), intent(inout) :: res
  integer(kind=c_int) :: return_code
  type(c_ptr), dimension(:), allocatable :: p
  integer(kind=c_int) :: num_optics
  integer :: i
  num_optics = size(optics)
  allocate(p(num_optics))
  do i = 1, num_optics
    p(i) = optics(i)%optics
  enddo
  return_code = c_add_optics(p, num_optics, res%optics)
  deallocate(p)
end function f_add_optics


function f_create_solar_flux(solar_flux, grid, path) &
  result(return_code)
  type(SolarFlux_t), intent(inout) :: solar_flux !< Solar flux object.
  type(Grid_t), intent(in) :: grid !< Spectral grid.
  character(kind=c_char, len=*), intent(in) :: path !< Solar flux csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(path, buf)
  solar_flux%solar_flux = c_null_ptr
  return_code = malloc_struct(solar_flux%solar_flux, solar_flux_struct)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = c_create_solar_flux(solar_flux%solar_flux, grid%grid, buf)
  deallocate(buf)
end function f_create_solar_flux


function f_destroy_solar_flux(solar_flux) &
  result(return_code)
  type(SolarFlux_t), intent(inout) :: solar_flux !< Solar flux object.
  integer(kind=c_int) :: return_code
  return_code = c_destroy_solar_flux(solar_flux%solar_flux)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = free_struct(solar_flux%solar_flux)
end function f_destroy_solar_flux


function f_solar_flux_properties(solar, flux) &
  result(return_code)
  type(SolarFlux_t), intent(in) :: solar !< Solar flux.
  real(kind=fp), dimension(:), intent(inout) :: flux !< Flux.
  integer(kind=c_int) :: return_code
  return_code = c_solar_flux_properties(solar%solar_flux, flux)
end function f_solar_flux_properties


function f_create_gas_optics(ml, num_levels, grid, device, hitran_path, h2o_ctm_dir, &
                             o3_ctm_dir, wcutoff, optical_depth_method) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: num_levels !< Number of atmospheric levels.
  type(Grid_t), intent(in) :: grid !< Spectral grid.
  type(Device_t), intent(in) :: device !< Device.
  character(kind=c_char, len=*), intent(in) :: hitran_path !< Path to HITRAN database file.
  character(kind=c_char, len=*), intent(in), optional :: h2o_ctm_dir !< Path to water vapor continuum directory.
  character(kind=c_char, len=*), intent(in), optional :: o3_ctm_dir !< Path to ozone continuum directory.
  real(kind=c_double), intent(in), optional :: wcutoff !< Cutoff from line center [1/cm].
  integer(kind=c_int), intent(in), optional :: optical_depth_method !< Method used to calculate the optical depths.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: hitran
  character(kind=c_char, len=1), dimension(:), allocatable :: h2ob
  character(kind=c_char, len=1), dimension(:), allocatable :: o3b
  call append_null_char(hitran_path, hitran)
  if (present(h2o_ctm_dir)) then
    call append_null_char(h2o_ctm_dir, h2ob)
  else
    call append_null_char("none", h2ob)
  endif
  if (present(o3_ctm_dir)) then
    call append_null_char(o3_ctm_dir, o3b)
  else
    call append_null_char("none", o3b)
  endif
  ml%ml = c_null_ptr
  return_code = malloc_struct(ml%ml, gas_optics_struct)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = c_create_gas_optics(ml%ml, num_levels, grid%grid, device%device, &
                                    hitran, h2ob, o3b, wcutoff, optical_depth_method)
  deallocate(hitran)
  deallocate(h2ob)
  deallocate(o3b)
end function f_create_gas_optics


function f_destroy_gas_optics(ml) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int) :: return_code
  return_code = c_destroy_gas_optics(ml%ml)
  if (return_code .ne. grtcode_success) then
    return
  endif
  return_code = free_struct(ml%ml)
end function f_destroy_gas_optics


function f_add_molecule(ml, molecule_id, min_line_center, max_line_center) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: molecule_id !< Molecule id.
  real(kind=c_double), intent(in), optional :: min_line_center !< Lower bound [1/cm] for spectral line centers.
  real(kind=c_double), intent(in), optional :: max_line_center !< Upper bound [1/cm] for spectral line centers.
  integer(kind=c_int) :: return_code
  return_code = grt_add_molecule(ml%ml, molecule_id, min_line_center, max_line_center)
end function f_add_molecule


function f_set_molecule_ppmv(ml, molecule_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: molecule_id  !< Molecule id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_molecule_ppmv(ml%ml, molecule_id, ppmv)
end function f_set_molecule_ppmv


function f_add_cfc(ml, cfc_id, filepath) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: cfc_id !< CFC id.
  character(kind=c_char, len=*), intent(in) :: filepath !< Path to CFC cross section csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(filepath, buf)
  return_code = grt_add_cfc(ml%ml, cfc_id, buf)
  deallocate(buf)
end function f_add_cfc


function f_set_cfc_ppmv(ml, cfc_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: cfc_id !< CFC id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_cfc_ppmv(ml%ml, cfc_id, ppmv)
end function f_set_cfc_ppmv


function f_add_cia(ml, species1, species2, filepath) &
  result(return_code)
  type(MolecularLines_t), intent(inout) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(in) :: species1 !< Id of species.
  integer(kind=c_int), intent(in) :: species2 !< Id of species.
  character(kind=c_char, len=*), intent(in) :: filepath !< Path to cross section csv file.
  integer(kind=c_int) :: return_code
  character(kind=c_char, len=1), dimension(:), allocatable :: buf
  call append_null_char(filepath, buf)
  return_code = grt_add_cia(ml%ml, species1, species2, buf)
  deallocate(buf)
end function f_add_cia


function f_set_cia_ppmv(ml, cia_id, ppmv) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecularlines object.
  integer(kind=c_int), intent(in) :: cia_id !< CIA species id.
  real(kind=fp), dimension(:), intent(in) :: ppmv !< Abundance [ppmv] (level).
  integer(kind=c_int) :: return_code
  return_code = grt_set_cia_ppmv(ml%ml, cia_id, ppmv)
end function f_set_cia_ppmv


function f_calculate_optics(ml, pressure, temperature, optics) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  real(kind=fp), dimension(:), intent(in) :: pressure !< Pressure [mb] (level).
  real(kind=fp), dimension(:), intent(in) :: temperature !< Temperature [K] (level).
  type(Optics_t), intent(inout) :: optics !< Optics object.
  integer(kind=c_int) :: return_code
  return_code = grt_calculate_optical_depth(ml%ml, pressure, temperature, optics%optics)
end function f_calculate_optics


function f_num_molecules(ml, n) &
  result(return_code)
  type(MolecularLines_t), intent(in) :: ml !< Molecular lines object.
  integer(kind=c_int), intent(out) :: n !< Number of molecules.
  integer(kind=c_int) :: return_code
  return_code = grt_get_num_molecules(ml%ml, n)
end function f_num_molecules


function f_grt_errstr(code, buf) &
  result(return_code)
  integer(kind=c_int), intent(in) :: code !< Error code.
  character(kind=c_char, len=*), intent(inout) :: buf !< Buffer to hold error message.
  integer(kind=c_int) :: return_code
  integer :: i
  logical :: eos
  return_code = c_grt_errstr(code, buf, len(buf))
  eos = .false.
  do i = 1, len(buf)
    eos = eos .or. buf(i:i) .eq. c_null_char
    if (eos) then
      buf(i:i) = " "
    endif
  enddo
end function f_grt_errstr


function f_rayleigh_scattering(optics, pressure) &
  result(return_code)
  type(Optics_t), intent(inout) :: optics !< Optics object.
  real(kind=fp), dimension(:), intent(in) :: pressure !< Pressure [mb] (level).
  integer(kind=c_int) :: return_code
  return_code = c_rayleigh_scattering(optics%optics, pressure)
end function f_rayleigh_scattering


end module grtcode
