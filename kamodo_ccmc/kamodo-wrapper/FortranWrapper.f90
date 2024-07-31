module FortranWrapper
use, intrinsic :: iso_c_binding
use,intrinsic :: iso_fortran_env
implicit none
interface
      type(c_ptr) function constructor(configpath_ptr) bind(c, name='constructor_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            character(len=:),pointer :: configpath_ptr
      end function constructor


      subroutine destructor(model_ptr) bind(c, name='destructor_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
      end subroutine destructor


      subroutine call_method_with_defaults(model_ptr, method_name_ptr, out_array_ptr, out_arr_shape, out_n_dim) &
                                    bind(c, name='call_method_with_defaults_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),pointer :: method_name_ptr
            type(c_ptr), value :: out_array_ptr
            type(c_ptr), value :: out_arr_shape
            integer(c_int),value :: out_n_dim
      end subroutine call_method_with_defaults

      subroutine call_method_with_input(model_ptr, method_name_ptr, in_array_ptr, in_arr_shape, in_n_dim,&
                  out_array_ptr, out_arr_shape, out_n_dim) bind(c, name='call_method_with_input_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),pointer :: method_name_ptr
            type(c_ptr), value :: in_array_ptr
            type(c_ptr), value :: in_arr_shape
            integer(c_int),value :: in_n_dim
            type(c_ptr), value :: out_array_ptr
            type(c_ptr), value :: out_arr_shape
            integer(c_int),value :: out_n_dim
      end subroutine call_method_with_input


      subroutine print_method_metadata(model_ptr, method_name_ptr) bind (c,name='print_method_metadata_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),pointer :: method_name_ptr
      end subroutine print_method_metadata

      subroutine print_method_defaults(model_ptr, method_name_ptr) bind (c,name='print_method_defaults_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),pointer :: method_name_ptr
      end subroutine print_method_defaults

      subroutine print_all_methods(model_ptr) bind (c,name='print_all_methods_KamodoWrapper')
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
      end subroutine print_all_methods


end interface

contains
      function constructor_kamodowrapper(configpath) result(model_ptr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: model_ptr
            character(len=:),allocatable,target :: configpath
            character(len=:),pointer::cf_path_ptr

            configpath = configpath//c_null_char
            cf_path_ptr => configpath
            model_ptr = constructor(cf_path_ptr)

      end function constructor_kamodowrapper


      subroutine destructor_kamodowrapper(model_ptr)
            use,intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr

            call destructor(model_ptr)
      end subroutine destructor_kamodowrapper


      subroutine call_method_with_defaults_kamodowrapper(model_ptr, method_name, out_array, out_arr_shape)
            use iso_c_binding
            implicit none
            type(c_ptr), value  :: model_ptr
            character(len=:),allocatable,target :: method_name
            real(c_double), dimension (:), allocatable, target :: out_array
            integer(c_int), dimension (:), target :: out_arr_shape
            character(len=:),pointer::method_name_ptr
            integer(c_int), dimension(:), pointer:: dummy_shape_out

            dummy_shape_out => out_arr_shape
            method_name = method_name //c_null_char
            method_name_ptr => method_name
            call call_method_with_defaults(model_ptr,method_name_ptr,c_loc(out_array),c_loc(dummy_shape_out),size(dummy_shape_out))

      end subroutine call_method_with_defaults_kamodowrapper


      subroutine call_method_with_input_kamodowrapper(model_ptr, method_name, in_array, in_arr_shape, out_array, out_arr_shape)
            use iso_c_binding
            implicit none
            type(c_ptr), value  :: model_ptr
            character(len=:),allocatable,target :: method_name
            real(c_double), dimension (:), allocatable, target  :: in_array, out_array
            integer(c_int), dimension (:), target :: in_arr_shape, out_arr_shape
            character(len=:),pointer::method_name_ptr
            integer(c_int), dimension(:), pointer:: dummy_shape_in, dummy_shape_out

            dummy_shape_in => in_arr_shape
            dummy_shape_out => out_arr_shape
            method_name = method_name//c_null_char
            method_name_ptr => method_name
            call call_method_with_input(model_ptr,method_name_ptr,c_loc(in_array),c_loc(dummy_shape_in),size(dummy_shape_in),&
                                                                  c_loc(out_array),c_loc(dummy_shape_out),size(dummy_shape_out))

      end subroutine call_method_with_input_kamodowrapper



      subroutine print_method_metadata_kamodowrapper(model_ptr, method_name)
            use iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),allocatable,target :: method_name
            character(len=:),pointer::method_name_ptr
            method_name=method_name//c_null_char
            method_name_ptr=>method_name
            call print_method_metadata(model_ptr, method_name_ptr)
      end subroutine print_method_metadata_kamodowrapper

      subroutine print_method_defaults_kamodowrapper(model_ptr, method_name)
            use iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            character(len=:),allocatable,target :: method_name
            character(len=:),pointer::method_name_ptr
            method_name=method_name//c_null_char
            method_name_ptr=>method_name
            call print_method_defaults(model_ptr, method_name_ptr)
      end subroutine print_method_defaults_kamodowrapper


      subroutine print_all_methods_kamodowrapper(model_ptr)
            use iso_c_binding
            implicit none
            type(c_ptr), value :: model_ptr
            call print_all_methods(model_ptr)
      end subroutine print_all_methods_kamodowrapper


end module FortranWrapper
