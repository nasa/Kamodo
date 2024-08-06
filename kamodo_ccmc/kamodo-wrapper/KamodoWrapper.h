#ifndef KAMODOWRAPPER_H
#define KAMODOWRAPPER_H


#include <string>
#include <map>
#include <Python.h>



class KamodoWrapper {
public:
    explicit KamodoWrapper(const std::string& configPath);
    virtual ~KamodoWrapper();

    void createKamodoInstance();
    void callMethodWithDefaults(const std::string& methodName, void* output, int* output_shape, int output_n_dim);
    void callMethodWithInput(const std::string& methodName, void* input, int* input_shape, int input_n_dim, void* output, int* output_shape, int output_n_dim);
    void printMethodMetaData(const std::string& methodName);
    void printMethodDefaults(const std::string& methodName);
    void printAllMethods();


private:



    std::string moduleName;
    PyObject* KamodoModule;
    std::string className;
    PyObject * KamodoClass;
    std::string configPath;
    PyObject * params;


    PyObject * kamodo_instance;
    PyObject * kamodo_signatures;
    std::map <std::string, PyObject*> variable_names; // Use py::object to store the list or any object
    PyObject * np_array_constructor;
    PyObject * np_array_reshape;
    PyObject * np_module;

    std::map<std::string, PyObject*> methodDefaults;
    std::map<std::string, PyObject *> methodMeta;

    void initializePython();
    void loadConfig(const std::string& configPath);
    PyObject* convert_to_numpy_array_and_unflatten(void* flat_array_ptr, int* desired_arr_shape, int desired_n_dim);
    void flatten_and_convert_to_c_array(void* array_ptr, PyObject* numpy_array, int* arr_shape, int n_dim);
};

extern "C" KamodoWrapper* constructor_KamodoWrapper(char ** configpath_ptr);
extern "C" void destructor_KamodoWrapper(void* model_ptr);
extern "C" void call_method_with_defaults_KamodoWrapper(void* model_ptr, char** methodName, void* output_array_ptr, int* output_arr_shape, int output_n_dim);
extern "C" void call_method_with_input_KamodoWrapper(void* model_ptr, char** methodName, void* input_array_ptr, int* input_arr_shape, int input_n_dim, void* output_array_ptr, int* output_arr_shape, int output_n_dim);
extern "C" void print_method_metadata_KamodoWrapper(void* model_ptr, char** methodName);
extern "C" void print_method_defaults_KamodoWrapper(void* model_ptr, char** methodName);
extern "C" void print_all_methods_KamodoWrapper(void* model_ptr);
#endif // KAMODOWRAPPER_H
