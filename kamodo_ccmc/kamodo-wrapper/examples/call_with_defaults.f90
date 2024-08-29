PROGRAM KamodoWrapperEx
use iso_c_binding
use FortranWrapper
implicit none




type(c_ptr) ::model_ptr
character(len=:),allocatable::cf_path, method_name

real(c_double), dimension(:,:,:), allocatable :: output

real(c_double), dimension(:), allocatable:: flat_output


cf_path = 'kamodo.yaml'

model_ptr = constructor_kamodowrapper(cf_path)

method_name = 'rho_ijk'

allocate(output(4,4,3))
allocate(flat_output(48))

call call_method_with_defaults_kamodowrapper(model_ptr,method_name, flat_output,shape(output))

output = reshape(flat_output,shape(output))
print *,output

deallocate(output)
deallocate(flat_output)

call destructor_kamodowrapper(model_ptr)

END PROGRAM KamodoWrapperEx
