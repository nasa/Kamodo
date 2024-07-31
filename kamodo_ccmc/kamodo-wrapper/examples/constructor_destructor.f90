PROGRAM KamodoWrapperEx
use iso_c_binding
use FortranWrapper
implicit none

type(c_ptr) ::model_ptr
character(len=:),allocatable::cf_path

cf_path = 'kamodo.yaml'
model_ptr = constructor_kamodowrapper(cf_path)

call destructor_kamodowrapper(model_ptr)

END PROGRAM KamodoWrapperEx
