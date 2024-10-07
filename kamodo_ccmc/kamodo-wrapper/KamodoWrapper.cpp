#include "KamodoWrapper.h"
#include <iostream>
#include <yaml-cpp/yaml.h>


//deallocated
KamodoWrapper::KamodoWrapper(const std::string& configPath) {
    Py_Initialize();



    this->loadConfig(configPath);
    this->createKamodoInstance(); // Optionally create the instance upon construction

    PyObject* np = PyImport_Import(PyUnicode_DecodeFSDefault("numpy"));
    PyObject *np_dict = PyModule_GetDict(np);
    PyObject *np_array_constructor =  PyDict_GetItemString(np_dict, "array");
    PyObject *np_array_reshape = PyDict_GetItemString(np_dict,"reshape");
    this-> np_array_constructor = np_array_constructor;
    this-> np_array_reshape = np_array_reshape;
    this->np_module = np;
}

//dealloacated
KamodoWrapper::~KamodoWrapper() {
    //remember to deallocate eveyrthing

    for (auto const& var : this->methodDefaults) {

        Py_DECREF(var.second);
    }
    for (auto const& var : this->methodMeta) {

        Py_DECREF(var.second);
    }


    Py_DECREF(this->np_module);


    Py_DECREF(this->params);
    Py_DECREF(this->kamodo_instance);
    Py_DECREF(this->KamodoClass);
    Py_DECREF(this->KamodoModule);
    Py_Finalize();

    std::cout << "Python interpreter finalized." << std::endl;
}


//deallocated
void KamodoWrapper::loadConfig(const std::string& configPath) {


    std::cout << "Loading configuration from: " << configPath << std::endl;
    // Load and process the YAML configuration file
        YAML::Node config = YAML::LoadFile(configPath);

        // Dynamically add module paths from the YAML file to sys.path
        if (config["module_paths"]) {
            for (const auto& path : config["module_paths"]) {
                std::string modulePath = path.as<std::string>();
                PyRun_SimpleString("import sys");
                std::string added_module_path= std::string("sys.path.append(\"")+modulePath+std::string("\")");
                PyRun_SimpleString(added_module_path.c_str());
            }
        }

        // Extract module and class name
        std::string modulePath = config["models"]["mymodel"]["class"].as<std::string>();
        size_t lastDot = modulePath.rfind(".");
        moduleName = modulePath.substr(0, lastDot);
        className = modulePath.substr(lastDot + 1);

        this->moduleName = moduleName;
        this->className = className;

        PyObject *kamodo_module = PyImport_Import(PyUnicode_DecodeFSDefault(this->moduleName.c_str()));
        if (PyErr_Occurred()) {
    // Handle the exception
    // Optionally, you can log the error or perform some cleanup
            PyErr_Print();
        }
        PyObject *dict = PyModule_GetDict(kamodo_module);
        PyObject *kamodo_class =  PyDict_GetItemString(dict, this->className.c_str());
        this->KamodoModule = kamodo_module;
        this->KamodoClass = kamodo_class;

        this->params = PyDict_New();

        for (YAML::const_iterator it = config["models"]["mymodel"]["params"].begin(); it != config["models"]["mymodel"]["params"].end(); ++it) {
            PyDict_SetItem(this->params,PyUnicode_FromString(it->first.as<std::string>().c_str()),PyUnicode_FromString(it->second.as<std::string>().c_str()));
        }
        Py_DECREF(dict);
}

//dealloacated
void KamodoWrapper::createKamodoInstance() {


    PyObject *kamodo_module = PyImport_Import(PyUnicode_DecodeFSDefault(this->moduleName.c_str()));
    PyObject *dict = PyModule_GetDict(kamodo_module);
    PyObject *kamodo_model =  PyDict_GetItemString(dict, this->className.c_str());



    //this is instantiate an instance w empty first argument and dictionary for keyword arguments
    PyObject *pArgs = PyTuple_Pack(0);
    this -> kamodo_instance = PyObject_Call(this->KamodoClass, pArgs,this->params);



    PyObject *model_dict = nullptr;
    model_dict = PyObject_GenericGetDict(this->kamodo_instance,nullptr);

    PyObject* signatures = PyDict_GetItemString(model_dict, "signatures");
    this->kamodo_signatures = signatures;

    PyObject  *pKey, *pValue,*var_ptr;
    Py_ssize_t pos = 0;


    while (PyDict_Next(this->kamodo_signatures, &pos, &pKey, &pValue)) {
        std::string var =std::string(PyBytes_AsString(PyUnicode_AsASCIIString(pKey)));
        var_ptr = PyObject_GetAttrString(this->kamodo_instance,var.c_str());
        this->variable_names.insert({var,var_ptr});
    }

    for (auto const& var : this->variable_names) {

        PyObject* var_meta = PyObject_GetAttrString(var.second,"meta");
        this->methodMeta.insert({var.first,var_meta});
    }

    PyObject* get_defaults_module = PyImport_Import(PyUnicode_DecodeFSDefault("get_defaults"));
    PyObject *get_defaults_dict = PyModule_GetDict(get_defaults_module);
    PyObject *get_default_args_func =  PyDict_GetItemString(get_defaults_dict, "get_default_args");

    for (auto const& var : this->variable_names) {
        PyObject* default_arg = PyObject_CallFunctionObjArgs(get_default_args_func,var.second,nullptr);
        this->methodDefaults.insert({var.first,default_arg});
    };
    Py_DECREF(get_defaults_dict);
    Py_DECREF(get_defaults_module);
    Py_DECREF(pKey);
    Py_DECREF(pValue);
    Py_DECREF(var_ptr);
    Py_DECREF(pArgs);
}


//dealloacated
PyObject* KamodoWrapper::convert_to_numpy_array_and_unflatten(void* flat_array_ptr, int* desired_arr_shape, int desired_n_dim) {
    int total_shape = 1;
    for (int i = 0; i < desired_n_dim; i++) {
        total_shape *= desired_arr_shape[i];
    }

    PyObject *dummy_list = PyList_New(total_shape);
    for (int i = 0; i < total_shape; i++) {
        PyObject* item = PyFloat_FromDouble(((double*)flat_array_ptr)[i]);
        PyList_SetItem(dummy_list, i, item);
    }
    PyObject* pArgs_np =  PyTuple_New(1);
    PyTuple_SetItem(pArgs_np, 0, dummy_list);

    PyObject *flat_np_arr = PyObject_CallMethod(this->np_module,"array","(O)",dummy_list);


    PyObject* pArgs_np_reshape = PyTuple_New(2);
    PyTuple_SetItem(pArgs_np_reshape, 0,flat_np_arr);
    PyObject* pArgs_dim =  PyTuple_New(desired_n_dim);

    for (int i = 0; i < desired_n_dim; i++) {
        PyTuple_SetItem(pArgs_dim, i, PyLong_FromLong(desired_arr_shape[i]));
    }
    PyTuple_SetItem(pArgs_np_reshape, 1,pArgs_dim);

    PyObject* fortran_ordering_reshape_arg = Py_BuildValue("{s:s}","order","F");
    //equivalent to calling np.reshape(flat_np_arr, desired_arr_shape,order='F')
    PyObject* reshaped_np_arr = PyObject_Call(this->np_array_reshape,pArgs_np_reshape,fortran_ordering_reshape_arg);


    Py_DECREF(dummy_list);
    Py_DECREF(pArgs_np);
    //cannot free the flat array here, because numpy.reshape doesnt make deep copy sometimes, so it can do weird stuff
    //Py_DECREF(flat_np_arr);
    Py_DECREF(pArgs_np_reshape);
    Py_DECREF(fortran_ordering_reshape_arg);

    return reshaped_np_arr;
}

//deallocated
void KamodoWrapper::flatten_and_convert_to_c_array(void* array_ptr, PyObject* numpy_array, int* arr_shape, int n_dim) {
    int total_shape = 1;
    for (int i = 0; i < n_dim; i++) {
        total_shape *= arr_shape[i];
    }

    PyObject* fortran_ordering_reshape_arg = Py_BuildValue("{s:s}","order","F");
    PyObject* flatten_numpy_array = PyObject_CallMethod(numpy_array,"flatten","(s)","F");
    PyObject* flatten_list = PyObject_CallMethod(flatten_numpy_array,"tolist",nullptr);
    for (int i = 0; i < total_shape; i++) {
        ((double*)array_ptr)[i] = PyFloat_AsDouble(PyList_GetItem(flatten_list,i));
    }
    Py_DECREF(fortran_ordering_reshape_arg);
    Py_DECREF(flatten_numpy_array);
    Py_DECREF(flatten_list);
}



//dealloacated
void KamodoWrapper::callMethodWithInput(const std::string& methodName, void* input, int* input_shape, int input_n_dim, void* output, int* output_shape, int output_n_dim) {
    if (this->variable_names.find(methodName)!=this->variable_names.end()) {
        PyObject* pfunc = this->variable_names[methodName];
        PyObject* input_numpy_array = this->convert_to_numpy_array_and_unflatten(input, input_shape, input_n_dim);
        PyObject *presult =  PyObject_Call(pfunc,Py_BuildValue("(O)",input_numpy_array) ,nullptr);
        this->flatten_and_convert_to_c_array(output, presult, output_shape, output_n_dim);
        Py_DECREF(presult);
        Py_DECREF(input_numpy_array);
    }
    else {
        std::cout<<"this variable doesn't exist\n";
        return;
    }
}

//deallocated
void KamodoWrapper::callMethodWithDefaults(const std::string& methodName, void* output, int* output_shape, int output_n_dim) {
    if (this->variable_names.find(methodName)!=this->variable_names.end()) {
        PyObject* pfunc = this->variable_names[methodName];
        PyObject *presult =  PyObject_Call(pfunc,Py_BuildValue("()") ,nullptr);
        this->flatten_and_convert_to_c_array(output, presult, output_shape, output_n_dim);
        Py_DECREF(presult);
    }
    else {
        std::cout<<"this variable doesn't exist\n";
        return;
    }

}

void KamodoWrapper::printMethodMetaData(const std::string& methodName) {
    if (this->methodMeta.find(methodName)!=this->methodMeta.end()) {
        PyObject_Print(methodMeta[methodName],stdout,0);
    }
    else {
        std::cout<<"this variable doesn't exist\n";
        return;
    }

}


void KamodoWrapper::printMethodDefaults(const std::string& methodName) {
    if (this->methodDefaults.find(methodName)!=this->methodDefaults.end()) {
        PyObject_Print(methodDefaults[methodName],stdout,0);
    }
    else {
        std::cout<<"this variable doesn't exist\n";
        return;
    }

}

void KamodoWrapper::printAllMethods() {
    for (auto const& var : this->variable_names) {
        std::cout<<var.first<<std::endl;
    };

}




extern "C" KamodoWrapper* constructor_KamodoWrapper(char ** configpath_ptr) {
    KamodoWrapper* model_ptr = new KamodoWrapper(*configpath_ptr);
    return model_ptr;
}

extern "C" void destructor_KamodoWrapper(void* model_ptr) {
    delete (KamodoWrapper*)model_ptr;
}


extern "C" void call_method_with_defaults_KamodoWrapper(void* model_ptr, char** methodName, void* output_array_ptr, int* output_arr_shape, int output_n_dim) {
    ((KamodoWrapper*)model_ptr)->callMethodWithDefaults(std::string(*methodName), output_array_ptr, output_arr_shape, output_n_dim);
}

extern "C" void call_method_with_input_KamodoWrapper(void* model_ptr, char** methodName, void* input_array_ptr, int* input_arr_shape, int input_n_dim, void* output_array_ptr, int* output_arr_shape, int output_n_dim) {
    ((KamodoWrapper*)model_ptr)->callMethodWithInput(std::string(*methodName), input_array_ptr, input_arr_shape, input_n_dim, output_array_ptr, output_arr_shape, output_n_dim);
}

extern "C" void print_method_metadata_KamodoWrapper(void* model_ptr, char** methodName) {
    ((KamodoWrapper*)model_ptr)->printMethodMetaData(std::string(*methodName));
}

extern "C" void print_method_defaults_KamodoWrapper(void* model_ptr, char** methodName) {
    ((KamodoWrapper*)model_ptr)->printMethodDefaults(std::string(*methodName));
}

extern "C" void print_all_methods_KamodoWrapper(void* model_ptr) {
    ((KamodoWrapper*)model_ptr)->printAllMethods();
}



// Implement other methods as needed...
