PROGRAM KamodoWrapperEx
use iso_c_binding
use FortranWrapper
implicit none

type(c_ptr) ::model_ptr
character(len=:),allocatable::cf_path, method_name
real(c_double), dimension(:,:), allocatable :: input
real(c_double), dimension(:,:,:), allocatable :: output
real(c_double), dimension(:), allocatable :: flat_input, flat_output
integer :: l, j


cf_path = 'kamodo.yaml'

model_ptr = constructor_kamodowrapper(cf_path)


allocate( input(3,3) )
allocate(flat_input(9))
!populate the input
do l = 1, 3
    do j = 1, 3
            input(l, j) = l*3 + j
    end do
end do
flat_input =reshape(input, shape(flat_input))

allocate(output(4,9,3))
allocate(flat_output(108))

method_name = 'rho_ijk'
call call_method_with_input_kamodowrapper(model_ptr, method_name,flat_input,shape(input),flat_output,shape(output))

output = reshape(flat_output,shape(output))
print *,output

deallocate(input)
deallocate(flat_input)
deallocate(output)
deallocate(flat_output)

call destructor_kamodowrapper(model_ptr)

END PROGRAM KamodoWrapperEx
