PROGRAM KamodoWrapperEx
use iso_c_binding
use FortranWrapper
implicit none




type(c_ptr) ::model_ptr
character(len=:),allocatable::cf_path, method_name



cf_path = 'kamodo.yaml'



model_ptr = constructor_kamodowrapper(cf_path)

print *, "all methods are:"
call print_all_methods_kamodowrapper(model_ptr)



method_name= 'rho_ijk'

print *,  method_name // "  metadata is: "
call print_method_metadata_kamodowrapper(model_ptr, method_name)

print *," "

print *, method_name // " default args are: "
call print_method_defaults_kamodowrapper(model_ptr, method_name)


print *," "
call destructor_kamodowrapper(model_ptr)

END PROGRAM KamodoWrapperEx
