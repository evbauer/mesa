! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module eos_wdm_av_load
   ! HDF5 loader for the WDM_AV processed per-isotope tables.  Reads files through
   ! forum (hdf5io_t), following the idiom of kap/private/kap_aesopus.f90 and
   ! rates/private/suzuki_tables.f90.  See eos/eosWDM_AV_builder/FORMAT.md.
   !
   ! Note on string attributes: forum reads character attributes with a
   ! fixed-length HDF5 string type, so the processed files must store the string
   ! attributes (isotope, element, ...) as fixed-length ASCII, not as h5py's
   ! default variable-length strings.

   use const_def, only: dp
   use eos_wdm_av_def
   use forum_m, only: hdf5io_t, OPEN_FILE_RO

   implicit none
   private

   public :: eos_wdm_av_read_table
   public :: eos_wdm_av_load_isotope
   public :: eos_wdm_av_free_tables

contains

   subroutine eos_wdm_av_read_table(filename, tbl, ierr)
      ! read one processed WDM_AV table from an HDF5 file into tbl
      character(len=*), intent(in) :: filename
      type(eos_wdm_av_table), intent(inout) :: tbl
      integer, intent(out) :: ierr

      type(hdf5io_t) :: hi
      integer :: fmt_version, incl_rad, has_zbar_i, has_eta_i
      logical :: file_exists

      ierr = 0

      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(*,'(a)') 'eos_wdm_av_read_table: file not found: '//trim(filename)
         ierr = wdm_av_ierr_bad_format
         return
      end if

      hi = hdf5io_t(filename, OPEN_FILE_RO)

      ! format version gate
      if (.not. hi% attr_exists('format_version')) then
         write(*,'(a)') 'eos_wdm_av_read_table: missing format_version attribute in '//trim(filename)
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if
      call hi% read_attr('format_version', fmt_version)
      if (fmt_version /= wdm_av_format_version) then
         write(*,'(a,i0,a,i0)') 'eos_wdm_av_read_table: unsupported format_version ', fmt_version, &
            ', this reader understands version ', wdm_av_format_version
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if

      ! radiation must be stripped from the processed tables
      if (.not. hi% attr_exists('includes_radiation')) then
         write(*,'(a)') 'eos_wdm_av_read_table: missing includes_radiation attribute'
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if
      call hi% read_attr('includes_radiation', incl_rad)
      if (incl_rad /= 0) then
         write(*,'(a)') 'eos_wdm_av_read_table: includes_radiation must be 0 (builder strips radiation)'
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if

      ! required scalar attributes
      if (.not. all_attrs_present(hi)) then
         write(*,'(a)') 'eos_wdm_av_read_table: missing a required grid attribute'
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if

      call hi% read_attr('Z', tbl% Z)
      call hi% read_attr('A', tbl% A)
      call hi% read_attr('interp_order', tbl% order)
      call hi% read_attr('nxi', tbl% nxi)
      call hi% read_attr('ntau', tbl% ntau)
      call hi% read_attr('xi_min', tbl% xi_min)
      call hi% read_attr('dxi', tbl% dxi)
      call hi% read_attr('tau_min', tbl% tau_min)
      call hi% read_attr('dtau', tbl% dtau)

      if (tbl% order /= wdm_av_order_bicubic .and. tbl% order /= wdm_av_order_biquintic) then
         write(*,'(a,i0)') 'eos_wdm_av_read_table: interp_order must be 3 or 5, got ', tbl% order
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if

      if (tbl% nxi < 2 .or. tbl% ntau < 2) then
         write(*,'(a,2i6)') 'eos_wdm_av_read_table: grid too small (need >= 2 per axis): ', tbl% nxi, tbl% ntau
         ierr = wdm_av_ierr_bad_format
         call hi% final()
         return
      end if

      ! optional zero-point offsets
      if (hi% attr_exists('e_offset')) call hi% read_attr('e_offset', tbl% e_offset)
      if (hi% attr_exists('f_offset')) call hi% read_attr('f_offset', tbl% f_offset)

      ! String attributes (isotope, element) are intentionally NOT read here:
      ! forum reads character attributes with a fixed-length native-char mem type,
      ! and HDF5 has no conversion path from the string types h5py writes (and a
      ! failed forum read aborts via error stop, which cannot be caught).  The
      ! isotope identity is set by the caller (eos_wdm_av_load_isotope) from the
      ! file name, matching the lookup policy in FORMAT.md section 3.

      ! required datasets (all shaped (nxi, ntau) in Fortran)
      call load_dset(hi, 'f', tbl% f, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xi', tbl% fx, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_tau', tbl% ft, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xixi', tbl% fxx, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_tautau', tbl% ftt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xitau', tbl% fxt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xixitau', tbl% fxxt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xitautau', tbl% fxtt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      call load_dset(hi, 'f_xixitautau', tbl% fxxtt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99

      ! optional Zbar block
      tbl% has_Zbar = .false.
      if (hi% attr_exists('has_Zbar')) then
         call hi% read_attr('has_Zbar', has_zbar_i)
         tbl% has_Zbar = (has_zbar_i /= 0)
      else
         tbl% has_Zbar = hi% dset_exists('Zbar')
      end if
      if (tbl% has_Zbar) then
         call load_dset(hi, 'Zbar', tbl% Zbar, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'Zbar_xi', tbl% Zbar_x, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'Zbar_tau', tbl% Zbar_t, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'Zbar_xitau', tbl% Zbar_xt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      end if

      ! optional eta block
      tbl% has_eta = .false.
      if (hi% attr_exists('has_eta')) then
         call hi% read_attr('has_eta', has_eta_i)
         tbl% has_eta = (has_eta_i /= 0)
      else
         tbl% has_eta = hi% dset_exists('eta')
      end if
      if (tbl% has_eta) then
         call load_dset(hi, 'eta', tbl% eta, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'eta_xi', tbl% eta_x, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'eta_tau', tbl% eta_t, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
         call load_dset(hi, 'eta_xitau', tbl% eta_xt, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      end if

      ! optional validity mask
      if (hi% dset_exists('valid')) then
         call load_dset_i1(hi, 'valid', tbl% valid, tbl% nxi, tbl% ntau, ierr); if (ierr /= 0) goto 99
      end if

      call hi% final()
      tbl% loaded = .true.
      return

99    continue
      call hi% final()
      tbl% loaded = .false.
   end subroutine eos_wdm_av_read_table


   logical function all_attrs_present(hi) result(ok)
      ! true iff every required grid attribute is present (each attr_exists call
      ! is made unconditionally, into a local, to avoid short-circuit warnings on
      ! the impure query function)
      type(hdf5io_t), intent(inout) :: hi
      logical :: e(9)
      e(1) = hi% attr_exists('Z')
      e(2) = hi% attr_exists('A')
      e(3) = hi% attr_exists('interp_order')
      e(4) = hi% attr_exists('nxi')
      e(5) = hi% attr_exists('ntau')
      e(6) = hi% attr_exists('xi_min')
      e(7) = hi% attr_exists('dxi')
      e(8) = hi% attr_exists('tau_min')
      e(9) = hi% attr_exists('dtau')
      ok = all(e)
   end function all_attrs_present


   subroutine load_dset(hi, name, arr, nxi, ntau, ierr)
      ! allocate-and-read a real(dp) dataset, verifying its (nxi, ntau) shape
      type(hdf5io_t), intent(inout) :: hi
      character(len=*), intent(in) :: name
      real(dp), allocatable, intent(inout) :: arr(:,:)
      integer, intent(in) :: nxi, ntau
      integer, intent(out) :: ierr
      ierr = 0
      if (.not. hi% dset_exists(name)) then
         write(*,'(a)') 'eos_wdm_av_read_table: missing dataset '//trim(name)
         ierr = wdm_av_ierr_bad_format
         return
      end if
      if (allocated(arr)) deallocate(arr)
      call hi% alloc_read_dset(name, arr)
      if (size(arr,1) /= nxi .or. size(arr,2) /= ntau) then
         write(*,'(a,2i7,a,2i7)') 'eos_wdm_av_read_table: dataset '//trim(name)//' has shape', &
            size(arr,1), size(arr,2), ' but expected', nxi, ntau
         ierr = wdm_av_ierr_bad_format
         return
      end if
   end subroutine load_dset


   subroutine load_dset_i1(hi, name, arr, nxi, ntau, ierr)
      ! allocate-and-read an integer(1) (int8) dataset, verifying its (nxi, ntau)
      ! shape.  forum has no int8 reader, so read into a default-integer buffer
      ! (HDF5 converts) and narrow to integer(1).
      type(hdf5io_t), intent(inout) :: hi
      character(len=*), intent(in) :: name
      integer(1), allocatable, intent(inout) :: arr(:,:)
      integer, intent(in) :: nxi, ntau
      integer, intent(out) :: ierr
      integer, allocatable :: buf(:,:)
      ierr = 0
      if (.not. hi% dset_exists(name)) then
         ierr = wdm_av_ierr_bad_format
         return
      end if
      call hi% alloc_read_dset(name, buf)
      if (size(buf,1) /= nxi .or. size(buf,2) /= ntau) then
         write(*,'(a,2i7,a,2i7)') 'eos_wdm_av_read_table: dataset '//trim(name)//' has shape', &
            size(buf,1), size(buf,2), ' but expected', nxi, ntau
         ierr = wdm_av_ierr_bad_format
         return
      end if
      if (allocated(arr)) deallocate(arr)
      allocate(arr(nxi, ntau))
      arr = int(buf, kind=1)
   end subroutine load_dset_i1


   subroutine eos_wdm_av_load_isotope(table_dir, isotope, idx, ierr)
      ! lazily load the table for isotope into the module store, returning its
      ! index.  Thread-safe via double-checked locking under OMP critical.
      character(len=*), intent(in) :: table_dir
      character(len=*), intent(in) :: isotope
      integer, intent(out) :: idx
      integer, intent(out) :: ierr
      character(len=256) :: filename

      ierr = 0
      idx = eos_wdm_av_table_index(isotope)
      if (idx > 0) return

      filename = trim(table_dir)//'/eosWDM_AV_'//trim(isotope)//'.h5'

!$OMP CRITICAL(eosWDM_AV_load)
      idx = eos_wdm_av_table_index(isotope)
      if (idx == 0) then
         if (.not. allocated(wdm_av_tables)) then
            allocate(wdm_av_tables(max_wdm_av_tables))
            num_wdm_av_tables = 0
         end if
         if (num_wdm_av_tables >= max_wdm_av_tables) then
            write(*,'(a,i0,a)') 'eos_wdm_av_load_isotope: table store is full (max ', &
               max_wdm_av_tables, ')'
            ierr = wdm_av_ierr_bad_format
         else
            call eos_wdm_av_read_table(filename, wdm_av_tables(num_wdm_av_tables+1), ierr)
            if (ierr == 0) then
               num_wdm_av_tables = num_wdm_av_tables + 1
               idx = num_wdm_av_tables
               if (len_trim(wdm_av_tables(idx)% isotope) == 0) &
                  wdm_av_tables(idx)% isotope = trim(isotope)
            else
               idx = 0
            end if
         end if
      end if
!$OMP END CRITICAL(eosWDM_AV_load)
   end subroutine eos_wdm_av_load_isotope


   subroutine eos_wdm_av_free_tables()
      ! release the module-level table store
      if (allocated(wdm_av_tables)) deallocate(wdm_av_tables)
      num_wdm_av_tables = 0
   end subroutine eos_wdm_av_free_tables

end module eos_wdm_av_load
